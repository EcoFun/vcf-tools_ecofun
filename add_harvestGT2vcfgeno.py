#!/usr/bin/env python2
# -*- coding: utf-8 -*-

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

#~#### VARIABLES AND LIBRARIES
#~# 0.1) import libraries
import argparse as ap, sys, gzip, vcf 

#~# 0.2) get options from command line
psr = ap.ArgumentParser(description="Add genotypes obtain from Harvest to Simon martin Geno files.")
psr.add_argument('harvest_vcf', help='vcf obtain from Harvest analysis. It contains three individuals only: the reference genome, the outgroup, and another instance of teh reference genome for control.', nargs=1)
psr.add_argument('geno', help='Geno file in Simon Martin format', nargs=1)
psr.add_argument('new_sample', help='Name of the new sample (for output header)', nargs=1)
psr.add_argument('new_vcf', help='output vcf including Harvest genotype', nargs=1)

#~psr.add_argument('-l', '--lparam', help='dummy default integer [40].', nargs=1, type=int, default=[40])
psr.add_argument('-a', '--allwrite', help='Verbose mode', action='store_true')
psr.add_argument('-v', '--verbose', help='Verbose mode', action='store_true')

#~# 0.3) set options
args = psr.parse_args()

fvcf = args.harvest_vcf[0]
geno = args.geno[0]
nom = args.new_sample[0]
out = args.new_vcf[0]
vb = args.verbose
al = args.allwrite
#~print al
#~sys.exit()
####### DEBBUG VARIABLES
#~fvcf = "Lprolificans-Sapioref_harvest.sorted.vcf.gz"
#~geno = "Scedo_26ind_AllPos.sort.recod.BI.geno.gz"
#~out = "00_data/Scedo_26ind_AllPos.sort.recod.hap.har.geno.gz"
#######

#~#### FUNCTIONS
def next_harv(vcfr):
	lh = next(vcfr)
	sh = lh.CHROM.split("_")[-1]
	ph = lh.POS
	return lh, sh, ph

#~#### MAIN SCRIPT
with gzip.open(geno) as g:
	lg = g.readline()
	with open(fvcf) as v:
		vcfr = vcf.Reader(v)
		lh, sh, ph = next_harv(vcfr)
		with gzip.open(out, "w") as o:
			while lg:
				if lg[0]=="#":
					o.write(" ".join([lg.strip(), nom]) + "\n")	# write output header
					lg = g.readline()	# set a geno SNP
				else:
					vg = lg.strip().split("\t")
					# scafs from both files
					sg = vg[0].split("_")[-1]
					sh = lh.CHROM.split("_")[-1]
					# check scaffolds are similar
					if int(sg) > int(sh):
						print "Scaffold of genotype file bigger than scaffold harvest vcf file!"
						lh, sh, ph = next_harv(vcfr)	# set a new harvest snp and restart loop
					elif int(sg) < int(sh):
						if al: o.write("\t".join([lg.strip(), "N"]) + "\n")	# write output
						lg = g.readline()	# set a new geno SNP and restart loop
					# else check pos
					else:
						if vb: print "Same scaf: %s" % sg
						pg = int(vg[1]) ; ph = lh.POS
						if vb: print pg, ph
						
						# check pos are the same
						if pg < ph:
							if vb: print "Same scaff, pass geno SNP"
							if al: o.write("\t".join([lg.strip(), "N"]) + "\n")	# write output
							lg = g.readline()	# set a new geno SNP and restart loop
						elif pg > ph:
							if vb: print "Same scaff, pass Harvest SNP"
							lh, sh, ph = next_harv(vcfr)	# set a new harvest snp and restart loop
						# if same pos, let's look at genotypes!:
						else:
							# if all fine, check harvest genotype
							gts = [ lh.samples[x]['GT'] for x in range(len(lh.samples)) ]
							if vb: print gts, lh.FILTER, lh.ALT
							if vb: print "Same scaf and position!"
							
							# check if len(ALT)==1, FILTER == PASS, and GT.ref == gt.contriol.ref
							if len(lh.FILTER) == 0 and len(lh.ALT)==1 and gts[0]==gts[-1]:
								if vb: print "Genotype written"
								o.write("\t".join([lg.strip(), str(lh.ALT[0])]) + "\n")	# write output
							else:
								if vb: print "Bad GT conditions :("
								if al: o.write("\t".join([lg.strip(), "N"]) + "\n")	# write output with N
							# do not forget to update SNPs!
							lg = g.readline()	# set a new geno SNP and restart loop
							lh, sh, ph = next_harv(vcfr)	# set a new harvest snp and restart loop



