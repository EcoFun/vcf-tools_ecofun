#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import argparse as ap, sys, gzip

# 0.1) get options from commands lines
parser = ap.ArgumentParser(description="Replace vcf 'GT' field wit ambiguous Allele depth by '.'.")

parser.add_argument('inp_vcf', help='vcf to process. BEWARE: the vcf must be sorted by coordinates!', nargs=1)
parser.add_argument('out_vcf', help='Corrected vcf.', nargs=1)

parser.add_argument('-d', '--min-dp', help='Minimum depth of sequencing for robust genotypes [5]', default=[6], nargs=1, type=int)
parser.add_argument('-o', '--out_bgt', help='File listing bad genotypes [list_bad_genotypes.csv.gz].', nargs=1, type=str, default=["list_bad_genotypes.csv.gz"])
parser.add_argument('-t', '--thres', help='Proportion of reads with the GT allele needed to validate GT [0.85].', nargs=1, type=float, default=[0.85])

# 0.2) set options
args = parser.parse_args()
vcfFile = args.inp_vcf[0]
fout = args.out_vcf[0]
fbad = args.out_bgt[0]
thres = args.thres[0]
mdp = args.min_dp[0]

#### FUNCTIONS
def open_fil(fil):
	if fil[-2:]=="gz":
		f = gzip.open(fil)
	else:
		f= open(fil)
	return(f)

def consistent_GT(GT, AD, mdp=mdp, sAD=sAD, thres=thres):
	test = False
	if sAD >=mdp:
		seuil = sAD-1 if sAD <= 10 else thres*sAD
		if GT == "0":
			test = int(AD[0]) >= seuil
		elif GT == "1":
			test = int(AD[1]) >= seuil
		else:
			sys.exit("Loci with more than two allele in the vcf")
	return test

### SCRIPT
# vcfFile can be gziped
vcf_reader = open_fil(vcfFile)
conOut = gzip.open(fout, "w")
conBad = gzip.open(fbad, "w")
count = 0
notInterest = 0
while True:
	try:
		l = vcf_reader.readline()
	except:
		print "\n#########\nFinished\n#########"
		break
	if l[0] == "#":
		conOut.write(l)
	else:
		ele = l.strip().split("\t")
		if ele[4] ==".":
			conOut.write(l)
		else:
			FORMAT = ele[8].split(":")
			iGT = FORMAT.index("GT")
			iAD = FORMAT.index("AD")
			for i, call in enumerate(ele[9:]):
				call = call.split(":")
				GT = call[iGT]
				AD = [ int(x) for x in call[iAD].split(",") ]
				sAD = sum(AD)
				test = consistent_GT(GT, AD)
				print call, test
				break
			break
				#~if GT == "." or :
					#~continue
				#~else:
					
	#~for individual in record:
		#~if individual["GT"] != ".":
			#~AlleleDepth = individual["AD"]
			#~if individual["GT"] == "1" and AlleleDepth[0] > AlleleDepth[1]:
				#~count += 1
				#~conBad.write("\t".join([individual.sample, record.CHROM, str(record.POS)])+"\n")
			#~else:
				
				#~notInterest = 1
				#~continue
	#~if notInterest:
		#~notInterest = 0
		#~continue

print count

conOut.close()
conBad.close()
