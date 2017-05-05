#!/usr/bin/env python

import vcf
import sys

# vcfFile can be gziped
vcfFile = sys.argv[1]

vcf_reader = vcf.Reader(open(vcfFile))


outFile = open("vcfProblems.list", "w")
count = 0
notInterest = 0
while True:
    try:
        record = vcf_reader.next()
    except:
        print "Finished"
        break
    if (record.ALT[0] is None):
        continue
    for individual in record:
        if individual["GT"] != ".":
            alternativesNumber = individual["AD"]
            if individual["GT"] == "1" and alternativesNumber[0] > alternativesNumber[1]:
                count += 1
                outFile.write("\t".join([individual.sample, record.CHROM, str(record.POS)])+"\n")
            else:
                notInterest = 1
                continue
    if notInterest:
        notInterest = 0
        continue

print count

outFile.close()
