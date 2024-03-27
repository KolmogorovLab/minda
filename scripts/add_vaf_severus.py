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
    sample_id = record.samples.keys()[0]
    record.info['VAF'] = record.samples[sample_id]['VAF']
    vcf_out.write(record)
vcf_out.close()