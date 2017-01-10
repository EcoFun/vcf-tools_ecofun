#!/usr/bin/env python2

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

import sys, argparse, vcf, Levenshtein as lv, itertools, numpy as np

# 0.1) get options from commands lines
parser = argparse.ArgumentParser(description=
'''Compute Dxy from a vcf. Currently the program removes SNPs with missing
data or with DP < DEPTH. Only two pops are allowed (specified in fpops). 
WARNING: the program requires a sorted vcf! See vcf-sort tool if needed. 
WARNING: Works only with haploid vcf.
WARNING: because the length of the last window of a scaffold cannot be 
known from a vcf, the program currently does not output result of the last 
window.''',
formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-a', '--allwind', 
	help='Output Dxy from all windows, including those with no SNPs.',
	action='store_true')
parser.add_argument('-d', '--depth', 
	help=
'''Minimum depth of sequencing for a SNP to be included
in the result file [10]''',
	nargs=1, default=[10], type=int)
parser.add_argument('-s', '--size', 
	help='Windows size in which to compute Dxy in bp [2000]',
	nargs=1, default=[2000], type=int)
parser.add_argument('-S', '--step', 
	help='Sliding step of the windows in bp. [SIZE]',
	nargs=1, default=[0], type=int)

parser.add_argument('fpops', help=
'''File listing the individuals belonging to pop1 and
pop2. Each pop is described on a single line,
individuals separeted by semi-colons.''', nargs=1)
parser.add_argument('vcf', help='vcf file', nargs=1)
parser.add_argument('output', help='output bed file', nargs=1)

# 0.2) set options
args = parser.parse_args()
#~print args
size = int(args.size[0])
step = int(args.step[0])
mdepth = int(args.depth[0])
awin = args.allwind

fpops = args.fpops[0]
finp = args.vcf[0]
fout = args.output[0]

if step > size:
	print "WARNING: window sliding step bigger than window size"

# if no sliding window
if step <1: step = size
#~print [size, step, mdepth, awin, fpops, finp, fout]
#~sys.exit()

#### DEBBUG
#~awin = False
#~size = 2000
#~step = 2000
#~mdepth = 10
#~fpops = 'pops_apio-boyd_dehooggi.txt'
#~fvcf = open("Scedo_26ind_twisstDxy.sorted.vcf")
#~out = open("test_dxy.bed", 'w')

#~fpops = 'pops_apio-boyd_dehooggi.txt'
#~finp = "Scedo_26ind_twisstDxy.sorted.vcf"
#~fout = "test_dxy.bed"
#### END of DEBBUG

# 0.3) functions
def set_pops(fpops):
	with open(fpops) as fp:
		pops = []
		for i, l in enumerate(fp):
			p = l.strip().split(";")
			pops.append(p)
		if len(pops[1])> len(pops[0]): 
			popsf = [ pops[1], pops[0] ]
		else:
			popsf = pops
		for j in range(len(popsf)):
			print "pop %s includes %s" % (j, popsf[j])
	return(popsf[0], popsf[1])

def set_coord(awin, pos, size, ch, pch, sto):
	# set initialization flag
	if not awin:	#info from SNP pos
		un = divmod(pos, size)
		# get coord
		if un[1] != 0:
			sta = un[0] * size + 1
			sto = (un[0] + 1) * size
		else:
			sto = pos
			sta = pos - size + 1
	elif ch != pch:	# all windows, new scaffold
		sta, sto = [1 , size]
	else:	# all windows, same scaff
		sta, sto = sto + 1, sto + size
	pch = ch	# in any cases
	return(sta, sto, pch)

def reset(awin, pos, size, ch, pch, sto, p0, p1):
	sta, sto, pch = set_coord(awin, pos, size, ch, pch, sto)
	ali = { 'p0': set_snp_dict(p0) , 'p1': set_snp_dict(p1) }
	lgw = size
	return(sta, sto, pch, ali, lgw)

def set_snp_dict(pop):
	mydict = {}
	for i in pop:
		 mydict [i] = ['']
	return(mydict)

def update_snp(rcd, p0, p1, mdepth, ch, pos, lgw, snp):
	for samp in rcd.samples:
		#~samp=rcd.samples[0]
		ind = samp.sample
		# test if ind in p0 or p1 (using try blocks as much quicker than if blocks)
		try:
			i = (p0+p1).index(ind)	# test for try
		except ValueError:
			#~print "Individual '%s' not in my pops" % ind
			continue	# so I stop here for this individual only!
		# once the test passed, let's go!
		#~print "upddate snp:", ind
		gt, dp = samp['GT'], samp['DP']
		# check if filters OK, good to check filters first here because 
			# there are not many individuals, would be too slow otherwise
		if gt[0] == "." or dp < mdepth:
			print "[MISSING DATA] missing data for individual '%s' at SNP '%s:%s'" % (ind, ch, pos)
			lgw -= 1
			return(lgw, None)	# I stop here for the whole SNP
		# if so, update snp given good pop 
		# (try blocks much quicker than if blocks)
		try:
			snp['p0'][ind][0] += gt
		except KeyError:
			snp['p1'][ind][0] += gt
	return(lgw, snp)

#### SCRIPT
p0, p1 = set_pops(fpops)

with open(fout, 'w') as out:
	with open(finp) as fvcf:
		# I. initialize output file
		vcf_reader = vcf.Reader(fvcf)
		noms = vcf_reader.samples
		head = "\t".join(["chrom", "start", "end", "N_SNP", "lg_win", "Dxy"])
		out.write(head+"\n")
		
		# II. initialize variables
		lgw = size
		ali = { 'p0': set_snp_dict(p0) , 'p1': set_snp_dict(p1) }
		
		# III. analyze vcf SNP per SNP
		#~rcd = next(vcf_reader)
		for rcd in vcf_reader:
			ch = rcd.CHROM
			pos = rcd.POS
			#~print "#####\n", ch, pos
			
			# III.1. if needed, set the first value of pch
			try:
				pch
			except NameError:
				sta, sto, pch = set_coord(awin, pos, size, ch, None, 0)
			
			# III.2. if change chr or new window, compute Dxy, print and 
			if ch != pch or pos > sto:
				# write only if not last window of the scaf
				if ch == pch:
					# if nothing in ali (i.e. no SNP)
					if len(ali['p0'][p0[0]][0]) == 0:
						# write only if output all windows
						if awin:
							res = [pch, str(sta), str(sto), "0", str(lgw), "0"]
							out.write("\t".join(res) + "\n")
					# if some SNPs are present:
					else:
						ckeys = list(itertools.product(ali['p0'].keys(), ali['p1'].keys()))
						ldist = []
						for c in ckeys:
							ldist.append(float(
								lv.distance(ali['p0'][c[0]][0], ali['p1'][c[1]][0])))
						dxy = round(sum(ldist)/(len(ldist)*lgw),5)
						# then print results (only if not last window of a scaffold)
						res = [pch, str(sta), str(sto), str(len(ali['p0'][p0[0]][0])), str(lgw), str(dxy) ]
						out.write("\t".join(res) + "\n")
				# reset variables in any case while new window
				sta, sto, pch, ali, lgw = reset(awin, pos, size, ch, pch, sto, p0, p1)
			
			# III.3 process the SNP: update alignment
			snp = { 'p0': set_snp_dict(p0) , 'p1': set_snp_dict(p1) }
			lgw, snp = update_snp(rcd, p0, p1, mdepth, ch, pos, lgw, snp)
			# if some data missing, then pass the current SNP
			if snp == None:
				#~print "PASS SNP"
				continue
			
			# if no missing data update ali (there was no break if code comes up to here)
			for k1 in snp.keys():
				for k2 in snp[k1].keys():
					ali[k1][k2][0] += snp[k1][k2][0]

	# DO NOT write the last window!
	#~out.write("\t".join(res) + "\n")

print '''
#### 
JOBE DONE
####
'''
