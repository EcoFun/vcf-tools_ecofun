#!/bin/env python2

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

#~from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys, argparse, os, gzip
from ete3 import Tree
from Bio.Alphabet import IUPAC


# 0.1) get options from commands lines
parser = argparse.ArgumentParser(description='Perform different manipulations on a fasta file. If no optional argument is provided, just format alignment to uppercase.')


parser.add_argument('outgroup', help='Compulsory outgroup.', nargs=1)
parser.add_argument('input_vcf', help='vcf to process.', nargs=1)
parser.add_argument('output_nwk', help='output files of newick trees.', nargs=1)

parser.add_argument('-k', '--keepali', help='Overwrite alignment files.', action='store_true')
parser.add_argument('-l', '--minlen', help='Minimum number of SNPs in order to keep windows on chromosomes edges [40].', nargs=1, type=int, default=[40])
parser.add_argument('-L', '--winlen', help='Length of the windows in SNPs [50].', nargs=1, type=int, default=[50])
parser.add_argument('-m', '--maxlenbp', help='Maximum window length in bp [5000].', nargs=1, type=int, default=[5000])
parser.add_argument('-M', '--mmiss', help='Maximum number of missing data per window for a given individual. In this case, the window is dropped (no tree computed here at all) [10].', nargs=1, type=int, default=[10])
parser.add_argument('-o', '--out_coord', help='Output file containing coordinates of retained windows [by default "output_nwk" with ".coord.csv" as extension].', nargs=1, type=str, default=[".coord.csv"])
parser.add_argument('-p', '--prefali', help='Prefix for alignment file(s) ["."].', nargs=1, type=str, default=["."])
parser.add_argument('-T', '--ntrees', help='Number of trees to compute [all]', nargs=1, type=int, default=[-1])
parser.add_argument('-u', '--uncompress', help='Overwrite alignment files [true].', action='store_false')
parser.add_argument('-v', '--verbose', help='verbose mode [false].', action='store_true')

	# options to add
		# minimum number of individuals per species => wait for Simon to implementthat in twisst
			# in order to remove some individuals with missing data

# 0.2) set options
args = parser.parse_args()
#~print args

outg = args.outgroup[0]
vcf = args.input_vcf[0]
nwk = args.output_nwk[0]
nmax = args.ntrees[0]
mmiss = args.mmiss[0]
lg = args.winlen[0]
lmax = args.maxlenbp[0]
minsnp = args.minlen[0]
prefali = args.prefali[0]
kali = args.keepali
vb = args.verbose
outc = args.out_coord[0]
if outc == ".coord.csv":
	outc = ".".join(nwk.split(".")[:-1]) + outc
#~print outc, nwk ; sys.exit()
# add gz extension if not present
if nwk[-2:] != "gz":
	nwk = nwk + ".gz"

#### FUNCTIONS
def print_info(verb, scaf, sta, sto, lgi, scaf2, nbad, bad):
	print "[ %s ] WINDOW %s:%s-%s" % (verb, scaf, sta, sto)
	print "number of SNPs: " + str(lgi)
	print scaf, scaf2
	print nbad
	print bad

def initialize(keys, ele):
	#~print "---"
	lgi = 1 ; seq = dict.fromkeys(keys)
	nbad = [0] * len(keys) ; bad = [False] * len(keys)
	scaf, sta = ele[0:2]
	return lgi, seq, nbad, bad, scaf, sta

def get_tree(n, keys, scaf, sta, sto, seq, outg, fnw, fout, lgi, vb):
	# increment tree number
	n += 1
	# prepare ali
	desc = "%s:%s-%s" % (scaf, sta, sto)
	if kali:
		fnom = "%s/ali.%s.%s.%s-%s.fasta" % (prefali, n, scaf, sta, sto)
	else:
		fnom = "%s/ali.%s.temp.fasta" % (prefali, n)
	with open(fnom, "w") as o:
		for k in keys:
			record = SeqRecord(Seq(seq[k], IUPAC.ambiguous_dna),
				id=k, description=desc)
			if vb: print record.seq
			o.write(record.upper().format("fasta"))
	
	# compute and retrieve tree using seaview...
	cmd = "seaview -build_tree -distance observed -NJ -o - %s" % fnom
	print cmd
	tr = os.popen(cmd).read().strip()
	if tr == '':
		n -= 1 
		cmd = "rm %s" % fnom
		os.system(cmd)
		return n
	tr = tr.split("] ")[1]
	# root tree
	tr = Tree(tr)
	tr.set_outgroup(outg)
	# write tree in a gz file
	tr = tr.write(format=1)
	fnw.write(tr + "\n")
	fout.write("\t".join([ str(x) for x in [scaf, sta, sto, (int(sta)+int(sto))/2,
		int(sto)-int(sta), lgi ] ]) + "\n")
	return n

#### SCRIPT
mess = "NORMAL EXIT: The program has computed %s trees\n" % nmax
with open(vcf) as f:
	with gzip.open(nwk, 'wb') as fnw:
		with gzip.open(outc, 'wb') as fout:
			fout.write("scaffold\tstart\tend\tmid\tlength\tsites\n")
			lgi, n = [1, 0]
			for l in f:
				l = l.strip()
				# if SNP line
				if l[0] != "#":
					ele = l.split("\t")
					if lgi == 1:
						print "---"
						scaf, sta = ele[0:2]
						sto = sta # just to avoid problem if only one snp on scaf
					else:
						scaf2 = ele[0]
						# check snp on same chromosome
						if scaf != scaf2:
							if lgi > minsnp:
								print_info("ALI SCAF EDGE", scaf, sta, sto, lgi, scaf2, nbad, bad)
								n = get_tree(n, keys, scaf, sta, sto, seq, outg, fnw, fout, lgi, vb)
								if n == nmax:
									sys.stderr.write(mess)
									break
								lgi, seq, nbad, bad, scaf, sta = initialize(keys, ele)
								lgi += 1
							else:
								print_info("PASS: DIFF SCAF", scaf, sta, sto, lgi, scaf2, nbad, bad)
								lgi, seq, nbad, bad, scaf, sta = initialize(keys, ele)
							continue
						
						sto = ele[1]	# should be updated only after scafs have been checked
						# check window not bigger than s in bp
						lw = int(sto)-int(sta)
						if lw > lmax:
							if lgi > minsnp:
								print_info("ALI LONG WINDOW", scaf, sta, sto, lgi, scaf2, nbad, bad)
								n = get_tree(n, keys, scaf, sta, sto, seq, outg, fnw, fout, lgi, vb)
								if n == nmax:
									sys.stderr.write(mess)
									break
								lgi, seq, nbad, bad, scaf, sta = initialize(keys, ele)
							else:
								print_info("PASS: LONG WINDOW", scaf, sta, sto, lgi, scaf2, nbad, bad)
								lgi, seq, nbad, bad, scaf, sta = initialize(keys, ele)
							continue

					# record allele info
					al = ele[3:5]
					for i,e in enumerate(ele[9:]):
						try:
							seq[keys[i]] += al[int(e)]
						except ValueError:
							try:
								seq[keys[i]] += "-"
								nbad[i] += 1
							except TypeError:
								seq[keys[i]] = "-"
								nbad[i] += 1
						except TypeError:
								seq[keys[i]] = al[int(e)]
					
					# check if no more than mmiss SNPs are missing
					for j,nb in enumerate(nbad):
						if nb > mmiss:
							bad[j] = True
					if True in bad:
						print_info("PASS: MISSING DATA", scaf, sta, sto, lgi, scaf2, nbad, bad)
						lgi, seq, nbad, bad, scaf, sta  = initialize(keys, ele)
						
						continue

					# if enough snps analysed, then write the temp fasta file
					if lgi == lg:
						print_info("EXPECTED ALI", scaf, sta, sto, lgi, scaf2, nbad, bad)
						n = get_tree(n, keys, scaf, sta, sto, seq, outg, fnw, fout, lgi, vb)
						if n == nmax:
							sys.stderr.write(mess)
							break
						lgi, seq, nbad, bad, scaf, sta  = initialize(keys, ele)
						
						# check if we reached the specified number of trees
					else:
						# increment number of analysed lines and proceed acccordingly.
						lgi += 1
				
				# pass header
				elif l[0:2] == "##":
					continue
				
				# if column header line, create dict: 1 entry per indiv (vcf header)...
				elif l[0:6] == "#CHROM":
					ele = l.split("\t")
					keys = ele[9:]
					seq = dict.fromkeys(keys)
					# initialize list of bad indiv
					nbad = [0] * len(keys) ; bad = [False] * len(keys)
					#~print seq
				else:
					sys.exit("ERROR: problem vcf format line:\n%s" % l)
