#!/usr/bin/python2.7
import sys, gzip, re

vcf = sys.argv[1]
outp = open(sys.argv[2], "w")
stop = sys.argv[3]

i = 0
sw = 0
with gzip.open(vcf) as f:
	for l in f:
		#~if "182571" in l: break
		if l[0] == "#": continue
		vcell = l.split("\t")
		alt = vcell[4]
		if alt == ".": continue	# SNPs only
		# get metrics name
		if sw == 0:
			info = vcell[7]
			linfo = re.split(";", info)
			lmet = [ x.split('=')[0] for x in linfo ]
			outp.write("\t".join(["QUAL"] + lmet) + "\n")
			sw = 1
		# extract metrics
		info = vcell[7]
		linfo = re.split(";|=", info)
		metrics = []
		#~break
		for x in lmet:
			try:
				metrics.append(linfo[linfo.index(x) + 1])
			except ValueError:
				metrics.append("NA")
		qual = vcell[5]
		res = "\t".join([qual] + metrics)
		#~break
		outp.write(res + "\n")
		i += 1
		if i % 5000 == 0:
			print i
			if i == stop:
				break

outp.close()

