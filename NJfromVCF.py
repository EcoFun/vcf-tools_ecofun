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
from re import split

# 0.1) get options from commands lines
parser = argparse.ArgumentParser(description='Compute NJ from vcf file.')

parser.add_argument('input_vcf', help='vcf to process. BEWARE: the vcf must be sorted per coordinates!', nargs=1)
parser.add_argument('output_nwk', help='output files of newick trees.', nargs=1)

parser.add_argument('-b', '--bed', help='Facultative bed file specifying window coordinates. BEWARE: the bedfile must be sorted per coordinates!', default=None, nargs=1)
parser.add_argument('-d', '--min-dp', help='Minimum depth of sequencing for robust genotypes [6]', default=[6], nargs=1, type=int)
parser.add_argument('-i', '--individuals', help='Individuals to process, separated by commas.', nargs=1, default=None)
parser.add_argument('-k', '--keepali', help='Keep all alignement files alignment files. If false, a singme temp file is used to write fasta files [False].', action='store_true')
parser.add_argument('-l', '--minlen', help='Minimum number of SNPs in order to keep windows on chromosomes edges [40].', nargs=1, type=int, default=[40])
parser.add_argument('-L', '--winlen', help='Length of sliding windows in SNPs [50].', nargs=1, type=int, default=[50])
parser.add_argument('-M', '--mmiss', help='Maximum number of missing data per individual per window. M is a SNP number for sliding window and a percentage for predefined windows. If missing data > M, the window is dropped [10].', nargs=1, type=int, default=[10])
parser.add_argument('-o', '--out_coord', help='Output file containing coordinates of retained windows [by default "output_nwk" with the ".coord.csv.gz" extension].', nargs=1, type=str, default=[".coord.csv.gz"])
parser.add_argument('-O', '--outgroup', help='Facultative outgroup for tree rooting.', nargs=1, default=None)
parser.add_argument('-p', '--prefali', help='Prefix for alignment file(s) ["."].', nargs=1, type=str, default=["."])
parser.add_argument('-s', '--split-pattern', help='Regular expression specifying a pattern used to split chromosomes name in order to obtain chromosome numbers (expected to be the last element of the chromosome name) ["_|chr|Chr"].', nargs=1, type=str, default=["_|chr|Chr"])
parser.add_argument('-T', '--ntrees', help='Number of trees to compute [all]', nargs=1, type=int, default=[-1])
parser.add_argument('-w', '--wsize', help='Maximum window length in bp [5000].', nargs=1, type=int, default=[5000])
parser.add_argument('-v', '--verbose', help='verbose mode [false].', action='store_true')

	# options to add
		# minimum number of individuals per species => wait for Simon to implementthat in twisst
			# in order to remove some individuals with missing data

# 0.2) set options
args = parser.parse_args()
#~print args

def open_fil(fil):
	if fil[-2:]=="gz":
		f = gzip.open(fil)
	else:
		f= open(fil)
	return(f)

bed = open_fil(args.bed[0]) if args.bed else None
outg = args.outgroup[0] if args.outgroup else None
outc = args.out_coord[0] ; nwk = args.output_nwk[0]
outc = ".".join(nwk.split(".")[:-1]) + outc if outc == ".coord.csv.gz" else outc
indivs = args.individuals[0].split(",") if args.individuals else None
mmiss = args.mmiss[0]
spatt = args.split_pattern[0]
vcf = args.input_vcf[0]

mdp = args.min_dp[0]
nmax = args.ntrees[0]
lg = args.winlen[0] if not bed else -1
wsize = args.wsize[0]
minsnp = args.minlen[0]
prefali = args.prefali[0]
kali = args.keepali
vb = args.verbose

# add gz extension if not present
if nwk[-2:] != "gz":
	nwk = nwk + ".gz"

#### FUNCTIONS
def print_info(verb, wscaf, wsta, wsto, lgi, sscaf, nbad, bad, pspos, suppl=None):
	print "---"
	print "[ %s ]\nWINDOW: %s:%s-%s" % (verb, wscaf, wsta, wsto)
	print "Number of SNPs: " + str(lgi)
	if suppl: print suppl
	print "SNP position: %s" % pspos
	print "Window and SNP scaffold: %s %s" % (wscaf, sscaf)
	print nbad
	print bad

def initialize_win(keys, l, bed=bed, wsize=wsize):
	ele = l.split("\t")
	lgi = 1 ; seq = dict.fromkeys(keys)
	nbad = [0] * len(keys) ; bad = [False] * len(keys)
	if not bed:
		try:
			wscaf, wsta = ele[0:2]
		except ValueError:
			return None, None, None, None, None, None, None, None
		wsto = str(int(wsta)+wsize)
	else:
		lbed = bed.readline()
		while lbed[0] == "#": lbed = bed.readline()
		v = lbed.strip().split("\t")
		wscaf, wsta, wsto = v[0:3]
		wsize = int(wsto) - int(wsta) + 1
	return lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize

def get_tree(n, keys, wscaf, wsta, wsto, seq, outg, fnw, fout, lgi, vb):
	# increment tree number
	n += 1
	# prepare ali
	desc = "%s:%s-%s" % (wscaf, wsta, wsto)
	if kali:
		fnom = "%s/ali.%s.%s.%s-%s.fasta" % (prefali, n, wscaf, wsta, wsto)
	else:
		fnom = "%s/ali.%s.temp.fasta" % (prefali, n)
	with open(fnom, "w") as o:
		for k in keys:
			record = SeqRecord(Seq(seq[k], IUPAC.ambiguous_dna),
				id=k, description=desc)
			if vb: print record.format("fasta").strip()
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
	tr = Tree(tr.split("] ")[1])
	# root tree
	if outg:
		tr.set_outgroup(outg)
	# write tree in a gz file
	tr = tr.write(format=1)
	fnw.write(tr + "\n")
	fout.write("\t".join([ str(x) for x in [wscaf, wsta, wsto, (int(wsta)+int(wsto))/2,
		int(wsto)-int(wsta), lgi ] ]) + "\n")
	return n

def test_leave(n, nmax):
	test = n == nmax
	if test:
		mess = "NORMAL EXIT: The program has computed %s trees\n" % nmax
		sys.stderr.write(mess)
	return test

def check_lgi_lg(wscaf, wsta, wsto, lgi, sscaf, nbad, bad, pspos, n, keys, 
	seq, outg, fnw, fout, vb):
	# output results and write log
	print_info("PASSED: MAX SNP WINDOW", wscaf, wsta, wsto, lgi, sscaf, nbad,
		bad, pspos)
	n = get_tree(n, keys, wscaf, wsta, wsto, seq, outg, fnw, fout, lgi, vb)
	return n

def check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto, sscaf, nbad, bad,
	n, keys, seq, outg, fnw, fout, vb, mess1, mess2, pspos, bed=bed, mmiss=mmiss):
	
	if bed:
		Msnp = round(((100.-mmiss)/100.) * lgi)
		Msnp = minsnp if Msnp < minsnp else Msnp
		vlgi = [ lgi - x for x in nbad ]
		vtest = [ x >= Msnp for x in vlgi ]
	else:
		Msnp = minsnp
		vtest = [ lgi >= Msnp ]
	
	sp = "Min SNP Nber: %s" % Msnp
	if not False in vtest:
		print_info(mess1, wscaf, wsta, wsto, lgi, sscaf, nbad, bad, pspos, suppl=sp)
		n = get_tree(n, keys, wscaf, wsta, wsto, seq, outg, fnw, fout, lgi, vb)
	else:
		print_info(mess2, wscaf, wsta, wsto, lgi, sscaf, nbad, bad, pspos, suppl=sp)
	return n

def rec_allele(alleles, vgts, seq, keys, nbad):
	for i,e in enumerate(vgts):
		call = e.strip().split(":")
		try:
			dp = int(call[2])
		except IndexError:	# there is no dp field
			dp = None
		except ValueError:	# dp is "."
			dp = None
		if dp and dp >= mdp:
			gt = call[0]
		else:
			gt = "bad"
		try:
			seq[keys[i]] += alleles[int(gt)]
		except ValueError:
			try:
				seq[keys[i]] += "N"
				nbad[i] += 1
			except TypeError:
				seq[keys[i]] = "N"
				nbad[i] += 1
		except TypeError:
				seq[keys[i]] = alleles[int(gt)]
	return seq, nbad

def update_SNP(f, spos):
	l = f.readline()
	pspos = spos
	return l, pspos

def proc_header(l, indivs, f):
	ele = l.strip().split("\t")
	keys = ele[9:] if not indivs else indivs
	inds = [ ele.index(x) for x in keys ]
	# then load data from first SNP
	l = f.readline()
	pspos = 0
	return keys, inds, l, pspos

def main(l, indiv):
	n = 0
	while 1:
		l = l.strip()
		# 1) PASSED vcf comment header
		if l[0:2] == "##":
			l = f.readline()	# just go to next line
			continue
		# 2) if vcf column header line, get sample names
		elif l[0:6] == "#CHROM":
			keys, inds, l, pspos = proc_header(l, indivs, f)
			continue
		# 3) if SNP line
		elif l[0] != "#":
			ele = l.strip().split("\t")
			sscaf, spos = ele[:2]
			# 3.1) if first SNP, initialize_win variables
			try:
				lgi
			except NameError:
				lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(keys, l)
			
			# 3.2) if we have the good number of SNPs in the window
			if not bed and lgi == lg:
				n = check_lgi_lg(wscaf, wsta, wsto, lgi, sscaf, nbad, bad, 
					pspos, n, keys, seq, outg, fnw, fout, vb)
				# re-initialize variables **from current SNP**
				lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(keys, l)
				# exit program if enough tree
				if test_leave(n, nmax): break
				# no need to restart the loop as a new window has been
					# initialized from the new SNP
			
			# 3.3) if SNP and window on different chromosomes: 
				# (to be checked BEFORE position checking)
			elif sscaf != wscaf:
				# 3.3.1) for sliding windows
					# no need to restart the loop as a new window has been
					# initialized from the new SNP
				if not bed:
					# output results if enough SNPs and write log
					n = check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto, 
						sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb,
						"PASSED: SCAFFOLD EDGE WINDOW, ENOUGH SNPs",
						"DISCARDED: SCAFFOLD EDGE WINDOW, NOT ENOUGH SNPs", pspos)
					# re-initialize variables **from current SNP**
					lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(keys, l)
					if test_leave(n, nmax): break
				
				# 3.3.2) for predefined windows
					# - In any case we change window, restarts loop to properly
						# process SNP as the new window can be far away.
					# - In any case we change SNP, restarts loop as the SNP can
						# still be on another chr or outside the window
				else:
					scnb = int(split(spatt, sscaf)[-1])
					wcnb = int(split(spatt, wscaf)[-1])
					# 3.3.2.1) SNP scaffold < window scaffold: update SNP
					if scnb < wcnb:
						l, pspos = update_SNP(f, spos)
						if not l: break
					# 3.3.2.2) SNP scaffold > window scaffold:
					elif scnb > wcnb:
						# output results if enough SNPs and write log
						n = check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto, 
							sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb,
							"PASSED: ENOUGH SNPS IN PREDEFINED WINDOW",
							"DISCARDED: NOT ENOUGH SNPS IN PREDEFINED WINDOW", pspos)
						if test_leave(n, nmax): break
						lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(keys, l)
					else:
						sys.exit("Scaff different but problem with ordering chr")
					continue

			# 3.4) if SNP position > window stop:
				# (to be done AFTER scaf checking)
			elif int(spos) > int(wsto):
				# output results if enough SNPs and write log
				if not bed:
					n = check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto,
						sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb,
						"PASSED: MAX LENGTH WINDOW, ENOUGH SNPs",
						"DISCARDED: MAX LENGTH WINDOW, NOT ENOUGH SNPs", pspos)
					# re-initialize variables **from current SNP**
					lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(keys, l)
					if test_leave(n, nmax): break
				else:
					n = check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto,
						sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb,
						"PASSED: ENOUGH SNPS IN PREDEFINED WINDOW",
						"DISCARDED: NOT ENOUGH SNPS IN PREDEFINED WINDOW", pspos)
					lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(keys, l)
					if test_leave(n, nmax): break
					continue	# restart loop if predefined windows
			
			# 3.5) check SNP position < window start	(only for predefined windows)
			elif int(spos) < int(wsta):
				if not bed:
					sys.exit("BUG: SNP position < window start with sliding window!")
				l, pspos = update_SNP(f, spos)
				if not l: break
				continue
			
			# 3.6) If allright, record allele info
			vgts = [ ele[x] for x in inds ]
			seq, nbad = rec_allele(ele[3:5], vgts, seq, keys, nbad)
			
			# 3.7) check if enough data for all individual in the current window
				# i.e. no more than 'miss' missing SNPs
			if not bed:
				for j,nb in enumerate(nbad):
					if nb > mmiss: bad[j] = True
				if True in bad:
					print_info("DISCARDED: TOO MUCH MISSING DATA PER INDIVIDUAL", wscaf, wsta, wsto, lgi,
						sscaf, nbad, bad, pspos)
					# thus change SNP and window!
					l, pspos = update_SNP(f, spos)
					if not l: break
					lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(keys, l)
					continue
			
			# 3.8) increment number of analysed lines and proceed acccordingly.
			l, pspos = update_SNP(f, spos)
			if not l: break
			lgi += 1
		else:
			sys.exit("ERROR: problem vcf format line:\n%s" % l)
	return (lgi, minsnp, wscaf, wsta, wsto, lgi, sscaf, nbad, bad, n, keys,
		seq, outg, fnw, fout, vb, pspos, ele)

#### SCRIPT
f = open_fil(vcf)
with gzip.open(nwk, 'wb') as fnw:
	with gzip.open(outc, 'wb') as fout:
		fout.write("scaffold\tstart\tend\tmid\tlength\tsites\n")
		l = f.readline()	# read vcf very first line
		# I) process vcf
		lgi, minsnp, wscaf, wsta, wsto, lgi, sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb, pspos, ele = main(l, indivs)

		# II) check last windows if enough snps analysed, then write the temp fasta file
		if not bed:
			n = check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto, 
				sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb,
				"LAST WINDOW PASSED: ENOUGH SNPS IN SLIDING WINDOW",
				"LAST WINDOW DISCARDED: NOT ENOUGH SNPS IN SLIDING WINDOW", pspos)
		else:
			n = check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto, 
				sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb,
				"LAST WINDOW PASSED: ENOUGH SNPS IN PREDEFINED WINDOW",
				"LAST WINDOW DISCARDED: NOT ENOUGH SNPS IN PREDEFINED WINDOW", pspos)

f.close()
if bed: bed.close()
