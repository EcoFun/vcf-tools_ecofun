#!/usr/bin/env python2

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

import sys, argparse, vcf 

# 0.1) get options from commands lines
parser = argparse.ArgumentParser(description='Compute density of heterozygous SNPs from a vcf')
parser.add_argument('-s', '--size', 
	help='Windows size in which to compute the GC in bp [1000]',
	nargs=1, default=[2000], type=int)
parser.add_argument('-S', '--step', 
	help='Sliding step of the windows in bp. By default the sliding step is "size" [size]',
	nargs=1, default=[0], type=int)
parser.add_argument('-D', '--depth', 
	help='Minimum depth of sequencing for a SNP to be included in the result file',
	nargs=1, default=[15], type=int)
parser.add_argument('vcf', help='vcf file', nargs=1)
parser.add_argument('output', help='output file', nargs=1)

# 0.2) set options
args = parser.parse_args()

size = int(args.size[0])
step = int(args.step[0])
mdepth = int(args.depth[0])

finp = args.vcf[0]
fout = args.output[0]

if step > size:
	print "WARNING: window sliding step bigger than window size"

# if no sliding window
if step <1: step = size
print args

#### DEBBUG
#~size = 2000
#~step = 2000
#~mdepth = 15
#~fvcf = open("01_Scedo_Diplo_Control_BadRegions.sorted.vcf")
#~out = open("test.txt", 'w')

#~finp = "01_Scedo_Diplo_Control_BadRegions.sorted.vcf"
#~fout = "test.txt"
#### END of DEBBUG

# 0.3) functions
def reset(ch, pch, size, step, sta, sto, coln):
	# if new chr, restart sta sto from start and change pch
	if ch != pch:
		print ch
		pch = ch
		sta, sto = [1, size]
	# if new window, just update sta and sto by step
	else:
		sta, sto = [ x + step for x in [ sta, sto ] ]
	# always reset counter in any case
	ct = [0] * len(coln)
	return(sta, sto, pch, ct)

#### SCRIPT
sta, sto = [1, size]
with open(fout, 'w') as out:
	with open(finp) as fvcf:
		#~rcd = next(vcf_reader)
		# initialize output file
		vcf_reader = vcf.Reader(fvcf)
		noms = vcf_reader.samples
		coln = [ '{}_{}'.format(a,b) for a in noms for b in ["homo", "hetero"] ]
		head = "\t".join(["chrom", "start", "end"] + coln )
		out.write(head+"\n")
		
		# initialize counter
		ct = [0] * len(coln)

		# analyze vcf SNP per SNP
		for rcd in vcf_reader:
			ch = rcd.CHROM
			pos = rcd.POS
			
			# if needed, set the first value of pch
			try:
				pch
			except NameError:
				pch = ch
			
			# if change chr or new window, print results and re-initialize variables
			if ch != pch or pos > sto:
				# print
				res = [pch, str(sta), str(sto)] + [ str(x) for x in ct ]
				out.write("\t".join(res) + "\n")
				# reset variables
				sta, sto, pch, ct = reset(ch, pch, size, step, sta, sto, coln)

			# then process all individuals at that SNP
			for i, samp in enumerate(rcd.samples):
				# check alleles are present and depth is Ok
				gt = samp['GT']
				dp = samp['DP']
				if gt[0] != "." and dp > mdepth:
					# if so, then count
					if gt[0] == gt[-1]:
						ct [i*2] += 1
					else:
						ct [i*2+1] += 1

	# DO NOT FORGET to write the last window!
	out.write("\t".join(res) + "\n")










