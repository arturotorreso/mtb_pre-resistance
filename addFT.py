#!/usr/bin/env python

import pysam
import sys
import vcf

input_vcf = sys.argv[1]
output_vcf = sys.argv[2]

myvcf = pysam.VariantFile(input_vcf, "r")


# Add the FT field to header.
if 'FT' not in myvcf.header.formats:
    myvcf.header.formats.add("FT", "1", "Integer", "Whether a sample was a Pass(1) or fail (0) based on FORMAT values")

with open(output_vcf, 'w') as output:
    output.write(str(myvcf.header))
    for variant in myvcf:
        for sample in variant.samples:
            gt = variant.samples[sample]['GT'][0]
            if 'PASS' not in variant.filter:
                variant.samples[sample]['FT'] = 0
                continue
            elif variant.samples[sample]['GT'][0] is None or variant.samples[sample]['GT'][0] == '.':
                variant.samples[sample]['FT'] = 0
                continue

            elif (variant.samples[sample]['GT'][0] == variant.samples[sample]['GT'][1]):
                if 'GQ' in variant.samples[sample]:
                    gq = variant.samples[sample]['GQ']
                else:
                    gq = 0
                dp = variant.samples[sample]['DP']
                sp = variant.samples[sample]['SP']
                alt = variant.samples[sample]['AD'][gt]
                falt = variant.samples[sample]['ADF'][gt]
                ralt = variant.samples[sample]['ADR'][gt]

                if sp is None:
                    sp = 0

                try:
                    falt_freq = variant.samples[sample]['ADF'][gt] / sum(variant.samples[sample]['ADF'])
                except ZeroDivisionError:
                    falt_freq = 0
                try:
                    ralt_freq = variant.samples[sample]['ADR'][gt] / sum(variant.samples[sample]['ADR'])
                except ZeroDivisionError:
                    ralt_freq = 0
                af1 = variant.samples[sample]['AD'][gt] / sum(variant.samples[sample]['AD'])

                if gq >= 0 and dp >= 10 and sp <= 50 and af1 >= 0.75 and falt_freq >= 0 and ralt_freq >= 0 and alt >= 10 and falt >= 2 and ralt >= 2:
                    variant.samples[sample]['FT'] = 1
                    continue
                else:
                    variant.samples[sample]['FT'] = 0
                    continue
            elif (variant.samples[sample]['GT'][0] != variant.samples[sample]['GT'][1]):
                gt_0 = variant.samples[sample]['GT'][0]
                gt_1 = variant.samples[sample]['GT'][1]

                if gt_0 != 0:
                    if variant.samples[sample]['AD'][gt_0] > variant.samples[sample]['AD'][gt_1]:
                        gt_1 = gt_0

                elif (variant.samples[sample]['GT'][0] == variant.samples[sample]['GT'][1]):
                    if 'GQ' in variant.samples[sample]:
                        gq = variant.samples[sample]['GQ']
                    else:
                        gq = 0
                dp = variant.samples[sample]['DP']
                sp = variant.samples[sample]['SP']
                alt = variant.samples[sample]['AD'][gt_1]
                falt = variant.samples[sample]['ADF'][gt_1]
                ralt = variant.samples[sample]['ADR'][gt_1]
                if sp is None:
                    sp = 0

                try:
                    falt_freq = variant.samples[sample]['ADF'][gt_1] / sum(variant.samples[sample]['ADF'])
                except ZeroDivisionError:
                    falt_freq = 0
                try:
                    ralt_freq = variant.samples[sample]['ADR'][gt_1] / sum(variant.samples[sample]['ADR'])
                except ZeroDivisionError:
                    ralt_freq = 0
                af1 = variant.samples[sample]['AD'][gt_1] / sum(variant.samples[sample]['AD'])

                if dp >= 10 and sp <= 50 and af1 > 0.75 and falt_freq > 0 and ralt_freq > 0 and alt >= 10 and falt >= 2 and ralt >= 2:
                    variant.samples[sample]['FT'] = 1
                    continue
                else:
                    variant.samples[sample]['FT'] = 0
                    continue
            else:
                variant.samples[sample]['FT'] = 0
        output.write(str(variant))
