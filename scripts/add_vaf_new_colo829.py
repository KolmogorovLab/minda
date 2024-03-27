#!/usr/bin/python3

import sys
import pysam
import statistics

vcf_file = sys.argv[1]

out_vcf = vcf_file.replace('.vcf' , '_vaf.vcf')
vcf_in=pysam.VariantFile(vcf_file,"r")
vcf_in.header.info.add("VAF",1,"Float","variant_allele_frequency")
vcf_out = pysam.VariantFile(out_vcf, 'w', header=vcf_in.header)

for record in vcf_in:
    support = list(record.info.values()[4])
    support = [sample.split("|") for sample in support]
    sample_vafs = [float(sample[2]) for sample in support]
    vaf = statistics.median(sample_vafs)
    record.info['VAF'] = vaf
    vcf_out.write(record)
vcf_out.close()