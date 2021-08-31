#!/usr/bin/env python

import sys
import argparse
import csv
import gzip
import binascii
import re
import os
from collections import OrderedDict
from typing import Dict
from warnings import warn


def main(input_vcf, output_maf, p_index, ref, wg, ins, dels, filts, merge_file, snps_only, gen_coords, hets, fail_as_ref, fail_as_N, fill_with_N, het_freq):

    # Check for type of input.

    if merge_file:
        merge_maf = merge_alignments(merge_file)
        with open(output_maf, 'w') as writer:
            for sample in merge_maf.keys():
                writer.write('>{}\n{}\n'.format(sample, merge_maf[sample]))

    elif input_vcf == '-':
        with sys.stdin as reader:
            # Parse file

            dictreader = _parseVCF(reader)

            # Write out file
            _getMAF(dictreader, output_maf, p_index, ref, wg, ins, dels, filts, snps_only, gen_coords, hets, fail_as_ref, fail_as_N, fill_with_N, het_freq)

    elif is_gz_file(input_vcf):
        with gzip.open(input_vcf, 'rt') as reader:
            # Parse file
            dictreader = _parseVCF(reader)

            # Write out file
            _getMAF(dictreader, output_maf, p_index, ref, wg, ins, dels, filts, snps_only, gen_coords, hets, fail_as_ref, fail_as_N, fill_with_N, het_freq)
    else:
        with open(input_vcf, 'rt') as reader:
            # Parse file
            dictreader = _parseVCF(reader)

            # Write out file
            _getMAF(dictreader, output_maf, p_index, ref, wg, ins, dels, filts, snps_only, gen_coords, hets, fail_as_ref, fail_as_N, fill_with_N, het_freq)


def is_gz_file(filepath):
    '''
    Return TRUE if file is gz compressed
    '''
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def _parseVCF(readable):
    '''
    Parse a VCF file into a dictionary with header as keys
    '''
    while True:
        line = readable.readline()
        # Check if the line is part of the upper header. If it is not, save it as header
        if not line.startswith('##'):
            header = line[1:].strip()
            break
    # Determine dialect
    dialect = 'excel-tab'
    # Header is the key fieldname for the dictionary
    fieldnames = header.split('\t')
    # Read file
    dictreader = csv.DictReader(readable, dialect=dialect, fieldnames=fieldnames)

    return dictreader


def _getMAF(dictreader, output_maf, p_index, ref, wg, ins, dels, filts, snps_only, gen_coords, hets, fail_as_ref, fail_as_N, fill_with_N, het_freq):

    samples = {key: [] for key in dictreader.fieldnames[9:]}
    ref_seq_parsed = parse_ref(ref)
    if fill_with_N:
        ref_seq_parsed = list('N'*len(ref_seq_parsed))

    index = []

    # Check for the validity of the data
    if len(samples) == 1 and ins:
        warn('\n** WARNING: Insertions are taken into account, but there is only one sample in the VCF file. This may lead to samples with different lengths if they are done independently **\n\n')
    if ins:
        warn('\n** Insertions currently not supported. Remove -i **\n\n')
        exit()
    if gen_coords:
        # Start from 0 or from the starting coordinate, until the end of the genome or the end coordinate
        prev_pos = int(gen_coords.split(':')[1].split('-')[0]) - 1
        end = int(gen_coords.split(':')[1].split('-')[1])
        if end > len(ref_seq_parsed):
            end = len(ref_seq_parsed)
    else:
        prev_pos = 0
        end = len(ref_seq_parsed)

    samples_starts = {key: prev_pos for key in dictreader.fieldnames[9:]}
    with open(p_index, 'w') as writer:
        n = 1
        # Write the index with variation
        writer.write('Alignment_position\tReference_position\n')
        last_pos = 0
        for line in dictreader:

            line = _parseVar(line, samples)
            pos = int(line['POS'])
            last_pos = pos + len(line['REF']) if pos + len(line['REF']) > last_pos else last_pos

            if pos < prev_pos:
                # Bit tricky: so the last variant from previous file was an indel that extended the end of previous file.
                # So we extend the start position here so we dont repeat.
                warn('\n** WARNING: An indel that started in a coordinate previous to the start point extends beyond the starting coordinate\nThis has been accounted for. The starting coordinate will be the end of that INDEL. **\n\n')
                rep_seq = pos + len(line['REF'])
                extra_seq = rep_seq - prev_pos - 1
                samples_starts.update((x, y + extra_seq) for x, y in samples_starts.items())
                continue
            if not _is_indel(line):
                # Analyze SNP
                index = pos
                writer.write('{}\t{}\n'.format(n, index))
                n += 1

                for sample in list(samples):

                    if pos <= samples_starts[sample]:
                        continue
                    if any(allele == '.' for allele in line[sample]['GT'][0].split('/')):
                        ref_kmer = ref_seq_parsed[samples_starts[sample]:pos - 1]
                        samples[sample] += ref_kmer
                        alt = resolve_iupac(line, sample)
                        samples[sample] += alt
                        samples_starts[sample] = pos
                        continue
                    # if all(int(allele) == 0 for allele in line[sample]['GT'][0].split('/')):
                    #     continue

                    if _filterGT(line, filts, sample) and _filterVar(line, filts):

                        ref_kmer = ref_seq_parsed[samples_starts[sample]:pos - 1]
                        samples[sample] += ref_kmer
                        if not is_het(line, sample):
                            gt = int(line[sample]['GT'][0].split('/')[0])
                            # alt = [line['ALT'].split(',')[gt - 1]] if line['ALT'] != '.' or gt != 0 else line['REF']
                            alt = [line['ALT'].split(',')[gt - 1]] if gt != 0 else line['REF']
                            samples[sample] += alt
                            samples_starts[sample] = pos

                        else:
                            alt = resolve_het(line, sample, hets, het_freq)
                            samples[sample] += alt
                            samples_starts[sample] = pos
                    else:
                        if fail_as_N:
                            ref_kmer = ref_seq_parsed[samples_starts[sample]:pos - 1]
                            samples[sample] += ref_kmer
                            alt = 'N'
                            samples[sample] += alt
                            samples_starts[sample] = pos

                        elif not fail_as_ref:
                            ref_kmer = ref_seq_parsed[samples_starts[sample]:pos - 1]
                            samples[sample] += ref_kmer
                            alt = resolve_iupac(line, sample)
                            samples[sample] += alt
                            samples_starts[sample] = pos

            if _is_indel(line) and (ins or dels):
                alts = line['ALT'].split(',')
                if (any(_is_del(alt, line['REF']) for alt in alts) and dels) or (any(_is_ins(alt, line['REF']) for alt in alts) and ins):
                    index = pos
                    writer.write('{}\t{}\n'.format(n, index))
                    n += 1
                for sample in list(samples):
                    if any(allele == '.' for allele in line[sample]['GT'][0].split('/')):
                        continue
                    if all(int(allele) == 0 for allele in line[sample]['GT'][0].split('/')):
                        continue
                    if _is_del(alts[0], line['REF']) and not dels:
                        continue
                    elif _is_ins(alts[0], line['REF']) and not ins:
                        continue
                    # if _filterGT(line, filts, sample) and _filterVar(line, filts) and any(_is_del(alt, line['REF']) for alt in alts):
                    if _filterGT(line, filts, sample) and _filterVar(line, filts):
                        ref_kmer = ref_seq_parsed[samples_starts[sample]:pos - 1]
                        samples[sample] += ref_kmer
                        if not is_het(line, sample):
                            gt = int(line[sample]['GT'][0].split('/')[0])
                            alt = alts[gt - 1]
                            alt = list(alt + '-'*(len(line['REF'])-len(alt)))
                            len_next = len(alt)
                            if pos <= samples_starts[sample]:
                                overlap = samples_starts[sample] - pos
                                alt = alt[overlap+1:]

                            samples[sample] += alt
                            samples_starts[sample] = pos + len_next - 1 if pos + len_next - 1 > samples_starts[sample] else samples_starts[sample]
                        else:
                            alt = resolve_het(line, sample, hets, het_freq)
                            len_next = len(alt)
                            if pos <= samples_starts[sample]:
                                overlap = samples_starts[sample] - pos
                                alt = alt[overlap+1:]

                            samples[sample] += alt
                            samples_starts[sample] = pos + len_next - 1 if pos + len_next - 1 > samples_starts[sample] else samples_starts[sample]
                    # elif any(_is_del(alt, line['REF']) for alt in alts):
                    else:
                        if not fail_as_ref:
                            ref_kmer = ref_seq_parsed[samples_starts[sample]:pos - 1]
                            samples[sample] += ref_kmer

                            alt = resolve_iupac(line, sample)
                            len_next = len(alt)
                            if pos <= samples_starts[sample]:
                                overlap = samples_starts[sample] - pos
                                alt = alt[overlap+1:]

                            samples[sample] += alt
                            samples_starts[sample] = pos + len_next - 1 if pos + len_next - 1 > samples_starts[sample] else samples_starts[sample]

    if output_maf:
        # Output to file if present
        with open(output_maf, 'w') as writer:
            if last_pos > end:
                end += (last_pos - end) - 1
            for sample in list(samples):
                ref_kmer = ref_seq_parsed[samples_starts[sample]:end]
                samples[sample] += ref_kmer
                writer.write('>{}\n{}\n'.format(sample, ''.join(samples[sample])))
    else:
        # Output to stdout if file not present
        if last_pos > end:
            end += (last_pos - end) - 1
        for sample in list(samples):
            ref_kmer = ref_seq_parsed[samples_starts[sample]:end]
            samples[sample] += ref_kmer
            print('>{}\n{}\n'.format(sample, ''.join(samples[sample])))


def _parseVar(var, samples):
    var['INFO'] = _parseINFO(var['INFO'])
    var_fmt = _parseFMT(var, samples)
    var['FILTER'] = _parseFILTER(var['FILTER'])
    return (var_fmt)


def _parseFILTER(var_filter):
    filt_lst = var_filter.split(';')
    return (filt_lst)


def _parseINFO(var_info):
    info_lst = var_info.split(';')
    if 'INDEL' in info_lst:
        info_lst[info_lst.index('INDEL')] = 'Type=INDEL'
    info_dict = {k: v for k, v in (x.split('=') for x in info_lst)}
    info_dict = {k: v.split(',') for k, v in info_dict.items()}
    return (info_dict)


def _parseFMT(var, samples):
    fmt = var['FORMAT'].split(':')
    for sample in list(samples):
        sample_lst = var[sample].split(':')
        sample_dict = {fmt[i]: sample_lst[i].split(',') for i in range(len(fmt))}
        var[sample] = sample_dict
    return (var)


def _filterVar(var, filts):
    '''Takes into acount QUAL, FILTER and INFO columns
    Returns True if filters are passed'''
    if not filts:
        return (True)
    var_filts = [filt for filt in filts if 'INFO' in filt or 'FILTER' in filt or 'QUAL' in filt]
    var_filts = [[x.strip() for x in re.split('(<=|>=|<|>|==|!=)', var_str)] for var_str in var_filts]
    for filt in var_filts:
        column = filt[0]
        operator = filt[1]
        val = filt[2]
        if val.isdigit():
            val = float(val)
        if column == 'FILTER':
            if eval("any(x {} '{}' for x in {})".format(operator, val, var[column])):
                return (False)
        elif 'INFO' in column:
            col = column.split('/')[0]
            fld = column.split('/')[1]
            index = re.search(r"\[([A-Za-z0-9_]+)\]", fld)
            if index:
                fld = re.sub("[\(\[].*?[\)\]]", "", fld)
                index = int(index.group(1))
            else:
                index = 0
            if isinstance(val, (int, float)):
                if eval("{} {} {}".format(var[col][fld][index], operator, val)):
                    return (False)
            else:
                if eval("{} {} '{}'".format(var[col][fld][index], operator, val)):
                    return (False)
        elif 'QUAL' in column:
            col = column
            if isinstance(val, (int, float)):
                if eval("{} {} {}".format(var[col], operator, val)):
                    return (False)

    return True


def _filterGT(var, filts, sample):
    '''Takes into acount FORMAT column
    Returns True if filters are passed'''
    if not filts:
        return (True)
    var_filts = [filt for filt in filts if 'FORMAT' in filt or 'FMT' in filt]
    var_filts = [[x.strip() for x in re.split('(<=|>=|<|>|==|!=)', var_str)] for var_str in var_filts]
    for filt in var_filts:
        fld = filt[0].split('/')[1]
        operator = filt[1]
        val = filt[2]
        try:
            val = float(val)
        except ValueError:
            val = val
        # if val.isdigit():
        #     val = float(val)
        index = re.search(r"\[([A-Za-z0-9_]+)\]", fld)
        if index:
            fld = re.sub("[\(\[].*?[\)\]]", "", fld)
            index = int(index.group(1))
        else:
            index = 0
        if var[sample][fld][index] == '.':
            return (False)
        if isinstance(val, (int, float)):
            if eval("{} {} {}".format(var[sample][fld][index], operator, val)):
                return (False)
        else:
            if eval("{} {} '{}'".format(var[sample][fld][index], operator, val)):
                return (False)
    return True


def _is_snp(alt, ref):
    if len(alt) == len(ref):
        return (True)
    return (False)


def _is_del(alt, ref):
    if len(alt) < len(ref):
        return (True)
    return (False)


def _is_ins(alt, ref):
    if len(alt) > len(ref):
        return (True)
    return (False)


def _is_indel(var):
    ref = var['REF']
    alts = var['ALT'].split(',')
    if len(ref) == 1 and all(len(alt) == 1 for alt in alts):
        return(False)
    return (True)


def is_het(var, sample):
    if len(var[sample]['GT'][0].split('/')) == 1:
        return False
    if any(gt == '.' for gt in var[sample]['GT'][0].split('/')):
        return False
    if (int(var[sample]['GT'][0].split('/')[0]) == int(var[sample]['GT'][0].split('/')[1])):
        return False
    return True


def resolve_iupac(var, sample):

    if _is_indel(var):
        alt_n = min(var['ALT'].split(','))
        alt_n = list(alt_n + 'N'*(len(var['REF'])-len(alt_n)))
        return(alt_n)
    else:
        if not is_het(var, sample):
            gt = int(var[sample]['GT'][0].split('/')[0]) if var[sample]['GT'][0].split('/')[0] != '.' else 0
            amb_alt = ''.join(sorted([var['REF'], var['ALT'].split(',')[gt - 1]])).upper() if var['ALT'] != '.' else 'GATC'
            if 'N' in amb_alt:
                amb_alt = 'GATC'
            return(list(ambiguous_dna_values.keys())[list(ambiguous_dna_values.values()).index(amb_alt)])
        else:
            if 'AD' not in var[sample]:
                return('N')
            amb_alt = [var['ALT'].split(',')[int(gt) - 1] for gt in var[sample]['GT'][0].split('/') if int(var[sample]['AD'][int(gt)])/int(var[sample]['DP'][0]) > 0.25 and int(gt) != 0]
            amb_alt = ''.join(sorted(amb_alt + [var['REF']])).upper()
            amb_alt = ''.join(sorted(list(dict.fromkeys(amb_alt))))
            if 'N' in amb_alt:
                amb_alt = 'GATC'
            return(list(ambiguous_dna_values.keys())[list(ambiguous_dna_values.values()).index(amb_alt)])


def resolve_het(var, sample, hets, het_freq):
    gt_0 = int(var[sample]['GT'][0].split('/')[0])
    gt_1 = int(var[sample]['GT'][0].split('/')[1])

    if 'AD' not in var[sample] and (hets == 'iupac' or hets == 'max'):
        hets = 'N'

    # if hets == 'N':
    #     if _is_indel(var):
    #         alt_n = min(var['ALT'].split(','))
    #         alt_n = list(alt_n + 'N'*(len(var['REF'])-len(alt_n)))
    #         return(alt_n)
    #     else:
    #         return('N')

    if hets == 'N':
        gt_pass = [int(gt) for gt in var[sample]['GT'][0].split('/') if int(var[sample]['AD'][int(gt)])/int(var[sample]['DP'][0]) > het_freq]
        if len(gt_pass) == 0:
            if _is_indel(var):
                alt_n = min(var['ALT'].split(','))
                alt_n = list(alt_n + 'N'*(len(var['REF'])-len(alt_n)))
                return(alt_n)
            else:
                return('N')
        elif len(gt_pass) > 1:
            if _is_indel(var):
                amb_alt = [var['REF'] if gt == 0 else var['ALT'].split(',')[gt - 1] for gt in gt_pass]
                max_alt = min(amb_alt)
                alt_n = list(max_alt + 'N'*(len(var['REF'])-len(max_alt)))
                # amb_alt = [var['REF'] if gt == 0 else var['ALT'].split(',')[gt - 1] for gt in gt_pass]
                # max_alt = max(amb_alt)
                # return ('N'*len(max_alt))
                return (alt_n)
            return('N')
        else:
            if _is_indel(var):
                alt = var['ALT'].split(',')[gt_pass[0] - 1]
                return (var['REF'] if gt_pass[0] == 0 else list(alt + '-'*(len(var['REF'])-len(alt))))
            else:
                return (var['REF'] if gt_pass[0] == 0 else var['ALT'].split(',')[gt_pass[0] - 1])

    if hets == 'iupac':
        gt_pass = [int(gt) for gt in var[sample]['GT'][0].split('/') if int(var[sample]['AD'][int(gt)])/int(var[sample]['DP'][0]) > het_freq]
        if len(gt_pass) == 0:
            if _is_indel(var):
                alt_n = min(var['ALT'].split(','))
                alt_n = list(alt_n + 'N'*(len(var['REF'])-len(alt_n)))
                return(alt_n)
            else:
                return('N')
        elif len(gt_pass) > 1:
            if _is_indel(var):
                amb_alt = [var['REF'] if gt == 0 else var['ALT'].split(',')[gt - 1] for gt in gt_pass]
                max_alt = min(amb_alt)
                alt_n = list(max_alt + 'N'*(len(var['REF'])-len(max_alt)))
                # amb_alt = [var['REF'] if gt == 0 else var['ALT'].split(',')[gt - 1] for gt in gt_pass]
                # max_alt = max(amb_alt)
                # return ('N'*len(max_alt))
                return (alt_n)
            amb_alt = ''.join(sorted([var['REF'] if gt == 0 else var['ALT'].split(',')[gt - 1] for gt in gt_pass], key=str.lower)).upper()
            if 'N' in amb_alt:
                amb_alt = 'GATC'
            return(list(ambiguous_dna_values.keys())[list(ambiguous_dna_values.values()).index(amb_alt)])
        else:
            if _is_indel(var):
                alt = var['ALT'].split(',')[gt_pass[0] - 1]
                return (var['REF'] if gt_pass[0] == 0 else list(alt + '-'*(len(var['REF'])-len(alt))))
            else:
                return (var['REF'] if gt_pass[0] == 0 else var['ALT'].split(',')[gt_pass[0] - 1])
    elif hets == 'max':
        max_ad = var[sample]['AD'].index(max(var[sample]['AD']))
        return (var['REF'] if max_ad == 0 else var['ALT'].split(',')[max_ad - 1])
    elif int(hets) == 1:
        return (var['REF'] if gt_0 == 0 else var['ALT'].split(',')[gt_1 - 1][0])
    elif int(hets) == 2:
        return ([var['ALT'].split(',')[gt_1 - 1]][0])


def parse_ref(ref):
    ref_list = []
    with open(ref, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                continue
            else:
                ref_list += list(line.rstrip())
    return (ref_list)


def parse_msa(msa):
    msa_dict = {}
    with open(msa, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                sample = line[1:].strip()
                msa_dict[sample] = []
                continue
            else:
                msa_dict[sample] += list(line.rstrip())
    return(msa_dict)


def merge_alignments(merge_file):
    with open(merge_file) as merge_files:
        # Check first file and create the dictionary
        first_file = next(merge_files).strip()
        merge_maf = parse_fasta(first_file)
        # Now go through the rest of alignments and add them to the correspoding kwy in the merge_maf
        for file_name in merge_files:
            maf_dict = parse_fasta(file_name.strip())
            for sample in merge_maf.keys():
                merge_maf[sample] += maf_dict[sample]
        return (merge_maf)


def parse_fasta(filename: str, ordered: bool = False) -> Dict[str, str]:
    """
    Parses a text file of genome sequences in fasta format (not fastq) into a dictionary.
    Arguments:
      filename: str - The name of the file containing the genome info.
      ordered: bool - Set this to True if you want the result to be ordered.
    """
    result = OrderedDict() if ordered else {}
    last_name = None
    with open(filename) as sequences:
        for line in sequences:
            if line.startswith('>'):
                last_name = line[1:].strip()
                result[last_name] = []
            else:
                result[last_name].append(line.strip())
    for name in result:
        result[name] = ''.join(result[name])
    return (result)


ambiguous_dna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "N": "GATC",
    }

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf_input', dest='input_vcf',
                        help='A VCF file. - for stdin')
    parser.add_argument('-o', '--output', dest='output_maf', default=False,
                        help='Output multialignment file (.maf)')
    parser.add_argument('-p', '--index', dest='p_index', default=os.devnull,
                        help='Index of reference and maf coordinates. Default: no file')
    parser.add_argument('-r', '--ref', dest='ref',
                        help='Reference genome used for the alignment. Not used unless -w/--whole.')
    parser.add_argument('-w', '--whole', dest='wg', action="store_true",
                        help='Output non variant sites. Needs the -r/--ref option')
    parser.add_argument('-d', '--dels', dest='dels', action="store_true", default=False,
                        help='''Include deletions in the output (noted as "-").
                        The output will have the same length as the initial reference genome.
                        ''')
    parser.add_argument('-i', '--ins', dest='ins', action="store_true", default=False,
                        help='''Include insertions in the output.
                        Output sequences will not have the same length as reference.
                        Gaps will be added into the samples without the insertions.
                        Only recommended with multi-sample VCFs.
                        In the future I will add a two-step process to add the insertions into single samples.
                        ''')
    parser.add_argument('-s', '--snps-only', dest='snps_only', action="store_true",
                        help='''Include only variant sites in the output alignment''')
    parser.add_argument('-c', '--coords', dest='gen_coords', default=False,
                        help='''
                        Coordinates in case that the input VCF doesnt include
                        the entire genome. Very relevant if only part of the genome
                        is given and the -w option is used.
                        Format: -c chromosome:start-end
                        E.g: bcftools view -c chromosome:start-end <vcf> | ./vcf2maf.py -v - -c chromosome:start-end
                        Default: uses all reference genome.
                        ''')
    parser.add_argument('-f', '--filt', dest='filts', action="append",
                        help='''
                        Filtering thresolds, for QUAL, FILTER, INFO and FORMAT columns. All the filters are to EXCLUDE.
                        One argument per filter.
                        Operators allowed: < | > | <= | >= | == | != \n
                        Eg: -f "FORMAT/DP < 5" -f "INFO/GQ < 20" -f "FILTER == PASS"
                        WARNING: For fields like FORMAT/AD you will need to specify
                        the index. This program uses 0-based indexes!!!!
                        Eg: for INFO/AD=200,10, INFO/AD[0] equals 200.
                        Look at the VCF specs for more info.
                        I guess I will expand it to other columns soon.
                        ''')
    parser.add_argument('-H', '--resolve-hets', dest='hets', choices=["1", "2", "max", "iupac", "N"], default='iupac',
                        help='''
                        Resolve heterozygous calls.
                        Arguments allowed: 1, 2, max, iupac or N.
                        1 if you want to always take the first allele (reference if 0/1).
                        2 if you want to always take the second allele (alternate if 0/1).
                        max if you want to take the allele with higher depth.
                        iupac use iupac ambiguity codes or alt/ref if the proportion of any of them is higher than 0.75.
                        N if you want all Hets as Ns.
                        IUPAC and Max requieres the tag FORMAT/AD. If not present, it will use 'N'.
                        Eg: --resolve-hets iupac
                        Default: iupac
                        ''')
    parser.add_argument('-t', '--het_freq', dest='het_freq', default=0.25,
                        help='''
                        Allele frequency to consider a heterozygous call. Lower values
                        bias the alignment towards more ambigous calls (N or IUPAC), higher values
                        towards a consensus alignment (unless it's too high). Default: 0.25
                        ''')

    parser.add_argument('-R', '--fail_as_ref', dest='fail_as_ref', action='store_true', default=False,
                        help="""Resolve low quality variants.
                        If present, low quality variants will show as the reference base
                        Otherwise, the iupac ambiguity code is used.
                        """)
    parser.add_argument('-N', '--fail_as_N', dest='fail_as_N', action='store_true', default=False,
                        help="""Resolve low quality variants.
                        If present, low quality variants will show as Ns.
                        Otherwise, the iupac ambiguity code is used.
                        """)
    parser.add_argument('-n', '--fill_with_N', dest='fill_with_N', action='store_true', default=False,
                        help="""If present, coordinates not present in the VCF file will be filled
                        with Ns instead of the reference genome.
                        """)
    parser.add_argument('-m', '--merge', dest='merge_file', default=False,
                        help='''File (absolute path please) with MAF files to merge, ordered by coordinate.
                        Overlaps between sequences is not allowed at the moment, so use for example 1-80 / 81-200.
                        Overlaps should be allowed and encouraged. Working on it at the moment.
                        Prefix files are assumed to have the same name as the .maf file,
                        with .prefix extension instead of .maf.
                        Use -o for output.
                        ''')

    args = parser.parse_args()

    if not len(sys.argv) > 1:
        parser.print_help()
        exit(1)

    if args.wg and args.snps_only:
        warn('\n** ERROR: You cant do a whole genome alignment and SNPs only! Come on, think about it! **\n\n')
        parser.print_help()
        exit(1)

    if args.ref is None:
            warn('\n** ERROR: reference fasta file required if -w/--whole flag present **\n\n')
            parser.print_help()
            exit(1)

    main(args.input_vcf,
         args.output_maf,
         args.p_index,
         args.ref,
         args.wg,
         args.ins,
         args.dels,
         args.filts,
         args.merge_file,
         args.snps_only,
         args.gen_coords,
         args.hets,
         args.fail_as_ref,
         args.fail_as_N,
         args.fill_with_N,
         float(args.het_freq)
         )
