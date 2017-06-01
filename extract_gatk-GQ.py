#!/usr/bin/python2.7
import sys, gzip, re

vcf = sys.argv[1]
outp = open(sys.argv[2], "w")
stop = sys.argv[3]

#~print vcf, outp
#~sys.exit()

i=0
outp.write("GQ" + "\n")
with gzip.open(vcf) as f:
	for l in f:
		if l[0] == "#": continue
		vcell = l.strip().split("\t")
		if "GQ" not in vcell[8]: continue	# SNPs only
		if i==0:
			igq = vcell[8].split(":").index("GQ")
			i += 1
		else:
			for c in vcell[9:]:
				if c == ".": continue
				res = 'NA' if c.split(":")[igq] == '.' else c.split(":")[igq]
				outp.write(res+"\n")
				#~print c.split(":")[igq]
				i += 1
			if i % 25000 == 0:
				print i
			if i >= stop:
				break

outp.close()
