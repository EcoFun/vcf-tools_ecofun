#!/usr/bin/env python2


# ~from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys, argparse as ap, os, gzip
from ete3 import Tree
from Bio.Alphabet import IUPAC
from re import split


__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

# 0.1) get options from commands lines
parser = ap.ArgumentParser(description='Compute NJ from vcf file.')

parser.add_argument('input_vcf', help='vcf to process. BEWARE: the vcf must be sorted by coordinates!', nargs=1)
parser.add_argument('output_nwk', help='output files of newick trees.', nargs=1)
parser.add_argument('-b', '--bed',
                    help='Facultative bed file specifying window coordinates. If the header has a header, it should start with a "#". BEWARE: the bedfile must be sorted per coordinates!',
                    default=None, nargs=1)
parser.add_argument('-d', '--min-dp',
                    help='Minimum depth of sequencing for robust genotypes [6]', default=[6], nargs=1, type=int)
parser.add_argument('-i', '--individuals',
                    help='Individuals to process, separated by commas.', nargs=1, default=None)
parser.add_argument('-k', '--keepali',
                    help='Keep all alignement files alignment files. If false, a temp file is used to write fasta files [False].', action='store_true')
parser.add_argument('-l', '--minlen',
                    help='Minimum number of good vcf entries to keep windows on chromosomes edges [40].', nargs=1, type=int, default=[40])
parser.add_argument('-L', '--winlen',
                    help='Length of sliding windows in SNPs [50]. Disabled if bed provided.', nargs=1, type=int, default=[50])
parser.add_argument('-M', '--mmiss',
                    help='Integer. Maximum number of missing data per individual per window. M is the number of discarded vcf entries while using sliding windows and a percentage while using predefined windows. If missing data > M, the window is dropped [10].', nargs=1, type=int, default=[10])
parser.add_argument('-N', '--notree',
                    help='Do not get tree (allow to print out fasta alignments using option -k) [False].', action='store_true')
parser.add_argument('-o', '--out_coord',
                    help='Output file containing coordinates of retained windows [by default "output_nwk" with the ".coord.csv.gz" extension].', nargs=1, type=str, default=[".coord.csv.gz"])
parser.add_argument('-O', '--outgroup',
                    help='Facultative outgroup for tree rooting.', nargs=1, default=None)
parser.add_argument('-p', '--prefali',
                    help='Prefix for alignment file(s) ["."].', nargs=1, type=str, default=["."])
parser.add_argument('-q', '--min-gtq',
                    help='Minimum "GQ" (genotype quality) required for robust genotypes [30]', default=[30], nargs=1, type=int)
parser.add_argument('-Q', '--min-qual',
                    help='Minimum "QUAL" (SNP quality) required for robust genotypes [30]', default=[30], nargs=1, type=int)
parser.add_argument('-s', '--split-pattern',
                    help='Regular expression specifying a pattern used to split chromosomes name in order to obtain chromosome numbers (expected to be the last element of the chromosome name) ["_|chr|Chr"].', nargs=1, type=str, default=["_|chr|Chr"])
parser.add_argument('-T', '--ntrees',
                    help='Number of trees to compute [all]', nargs=1, type=int, default=[-1])
parser.add_argument('-w', '--wsize',
                    help='Maximum window length in bp [5000].', nargs=1, type=int, default=[5000])
parser.add_argument('-v', '--verbose',
                    help='verbose mode [false].', action='store_true')

# options to add
# minimum number of individuals per species => wait for Simon to implementthat in twisst
# in order to remove some individuals with missing data

# 0.2) set options
args = parser.parse_args()
# ~print args


def open_fil(fil):
    if fil[-2:] == "gz":
        f = gzip.open(fil)
    else:
        f = open(fil)
    return(f)


bed = open_fil(args.bed[0]) if args.bed else None
outg = args.outgroup[0] if args.outgroup else None
outc = args.out_coord[0]
nwk = args.output_nwk[0]
outc = ".".join(nwk.split(".")[:-1]) + outc if outc == ".coord.csv.gz" else outc
indivs = args.individuals[0].split(",") if args.individuals else None
mmiss = args.mmiss[0]
spatt = args.split_pattern[0]
vcf = args.input_vcf[0]

mdp = args.min_dp[0]
mgtq = args.min_gtq[0]
mqual = args.min_qual[0]
nmax = args.ntrees[0]
lg = args.winlen[0] if not bed else -1
wsize = args.wsize[0]
minsnp = args.minlen[0]
prefali = args.prefali[0]
kali = args.keepali
notree = args.notree
vb = args.verbose

# add gz extension if not present
if nwk[-2:] != "gz":
    nwk = nwk + ".gz"

# FUNCTIONS


def print_info(verb, wscaf, wsta, wsto, lgi, sscaf, nbad, bad, pspos, vlgi=None, vtest=None, suppl=None):
    print "---"
    print "[ %s ]\nWINDOW: %s:%s-%s" % (verb, wscaf, wsta, wsto)
    print "Number of vcf entries processed: " + str(lgi)
    if suppl:
        print suppl
    print "Position of last vcf entry: %s" % pspos
    print "Window and vcf entry scaffold: %s %s" % (wscaf, sscaf)
    print "Number bad positions", nbad
    print "test maximum number of bad positions:", bad
    if vlgi:
        print "Number of valid vcf entries:", vlgi
    if vtest:
        print "Test minimum number of vcf entries:", vtest


def initialize_win(keys, l, bed=bed, wsize=wsize):
    ele = l.split("\t")
    lgi = 1
    seq = dict.fromkeys(keys)
    nbad = [0] * len(keys)
    bad = [False] * len(keys)
    if not bed:
        try:
            wscaf, wsta = ele[0:2]
        except ValueError:
            print """
####
VCF file totally processed
####"""
            return None, None, None, None, None, None, None, None
        wsto = str(int(wsta) + wsize)
    else:
        lbed = bed.readline()
        if not lbed:
            print """
####
Bed file totally processed
####"""
            return None, None, None, None, None, None, None, None
        while lbed[0] == "#":
            lbed = bed.readline()
        v = lbed.strip().split("\t")
        wscaf, wsta, wsto = v[0:3]
        wsize = int(wsto) - int(wsta) + 1
    return lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize


def get_tree(n, keys, wscaf, wsta, wsto, seq, outg, fnw, fout, lgi, vb, kali=kali, notree=notree):
	# increment tree number
	n += 1
	# prepare ali
	desc = "%s:%s-%s" % (wscaf, wsta, wsto)
	if kali:
		fnom = "%s/%s.%s-%s.fasta" % (prefali, wscaf, wsta, wsto)
	else:
		fnom = "%s/ali.temp.fasta" % (prefali)
	
	with open(fnom, "w") as o:
		for k in keys:
			record = SeqRecord(Seq(seq[k], IUPAC.ambiguous_dna),
				id=k, description=desc)
			if vb: print record.format("fasta").strip()
			o.write(record.upper().format("fasta"))

	if notree:
		return n

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
                 seq, outg, fnw, fout, vb, notree=notree):
    # output results and write log
    print_info("PASSED: MAX VCF ENTRY WINDOW", wscaf, wsta, wsto, lgi, sscaf, nbad,
               bad, pspos)
    n = get_tree(n, keys, wscaf, wsta, wsto, seq, outg, fnw, fout, lgi, vb)
    return n


def check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto, sscaf, nbad, bad,
                     n, keys, seq, outg, fnw, fout, vb, mess1, mess2, pspos, bed=bed,
                     mmiss=mmiss, notree=notree):

    if bed:
        Msnp = round(((100. - mmiss) / 100.) * lgi)
        Msnp = minsnp if Msnp < minsnp else Msnp
        vlgi = [lgi - x for x in nbad]
        vtest = [x >= Msnp for x in vlgi]
    else:
        Msnp = minsnp
        vlgi = lgi
        vtest = [vlgi >= Msnp]

    sp = "Min Nber of valid vcf entries: %s" % Msnp
    if not False in vtest:
        print_info(mess1, wscaf, wsta, wsto, lgi, sscaf, nbad,
                   bad, pspos, vlgi=vlgi, vtest=vtest, suppl=sp)
        n = get_tree(n, keys, wscaf, wsta, wsto, seq,
                     outg, fnw, fout, lgi, vb, vtest)
    else:
        print_info(mess2, wscaf, wsta, wsto, lgi, sscaf, nbad,
                   bad, pspos, vlgi=vlgi, vtest=vtest, suppl=sp)
    return n


def get_format(FORMAT, formID, iform):
    try:
        ind = FORMAT.index(formID)
        try:
            res = int(iform[ind])
        except ValueError:    # data at the field is "."
            res = "."
    except ValueError:    # FORMAT field not present
        res = None
    return res


def get_gt(iform, FORMAT, QUAL, mqual=mqual, mdp=mdp, mgtq=mgtq):
    dp = get_format(FORMAT, "DP", iform)
    gtq = get_format(FORMAT, "GQ", iform)
    if QUAL >= mqual and dp >= mdp and (gtq >= mgtq or not gtq):
        gt = iform[0]    # GT always the first field and always present.
    else:
        gt = "N"
    return gt


def rec_allele(vgts, seq, keys, alleles, nbad, FORMAT, QUAL,
               mqual=mqual, mdp=mdp, mgtq=mgtq):
    for i, e in enumerate(vgts):
        iform = e.strip().split(":")
        gt = get_gt(iform, FORMAT, QUAL, mqual, mdp, mgtq)
        try:
            # try to extend the sequence
            seq[keys[i]] += alleles[int(gt)]
        except ValueError:
            # GT is "N", so let's directly add 'N' to the sequence
            try:
                # try to extend the sequence with 'N'
                seq[keys[i]] += "N"
                nbad[i] += 1
            except TypeError:
                # the sequence is not initialized yet, let's do it
                seq[keys[i]] = "N"
                nbad[i] += 1
        except TypeError:
            # the sequence is not initialized yet, let's do it
            seq[keys[i]] = alleles[int(gt)]
    return seq, nbad


def update_SNP(f, spos):
    l = f.readline()
    pspos = spos
    return l, pspos


def proc_header(l, indivs, f):
    ele = l.strip().split("\t")
    keys = ele[9:] if not indivs else indivs
    inds = [ele.index(x) for x in keys]
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
            l = f.readline()    # just go to next line
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
                lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(
                    keys, l)
                if not lgi:
                    break

            # 3.2) if we have the good number of SNPs in the window
            if not bed and lgi == lg:
                n = check_lgi_lg(wscaf, wsta, wsto, lgi, sscaf, nbad, bad,
                                 pspos, n, keys, seq, outg, fnw, fout, vb)
                # re-initialize variables **from current SNP**
                lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(
                    keys, l)
                # exit program if enough tree
                if test_leave(n, nmax) or not lgi:
                    break
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
                                         "PASSED: SCAFFOLD EDGE WINDOW, ENOUGH VCF ENTRIES",
                                         "DISCARDED: SCAFFOLD EDGE WINDOW, NOT ENOUGH VCF ENTRIES", pspos)
                    # re-initialize variables **from current SNP**
                    lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(
                        keys, l)
                    if test_leave(n, nmax) or not lgi:
                        break

                # 3.3.2) for predefined windows
                # - In any case we change window, restarts loop to properly
                # process SNP as the new window can be far away.
                # - In any case we change SNP, restarts loop as the SNP can
                # still be on another chr or outside the window
                else:
                    print spatt
                    print sscaf
                    scnb = int(split(spatt, sscaf)[-1])
                    wcnb = int(split(spatt, wscaf)[-1])
                    # 3.3.2.1) SNP scaffold < window scaffold: update SNP
                    if scnb < wcnb:
                        l, pspos = update_SNP(f, spos)
                        if not l:
                            break
                    # 3.3.2.2) SNP scaffold > window scaffold:
                    elif scnb > wcnb:
                        # output results if enough SNPs and write log
                        n = check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto,
                                             sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb,
                                             "PASSED: ENOUGH VCF ENTRIES IN PREDEFINED WINDOW",
                                             "DISCARDED: NOT ENOUGH VCF ENTRIES IN PREDEFINED WINDOW", pspos)
                        if test_leave(n, nmax):
                            break
                        lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(
                            keys, l)
                        if not lgi:
                            break
                    else:
                        sys.exit(
                            "Scaff different but problem with ordering chr")
                    continue

            # 3.4) if SNP position > window stop:
                # (to be done AFTER scaf checking)
            elif int(spos) > int(wsto):
                # output results if enough SNPs and write log
                if not bed:
                    n = check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto,
                                         sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb,
                                         "PASSED: MAX LENGTH WINDOW, ENOUGH VCF ENTRIES",
                                         "DISCARDED: MAX LENGTH WINDOW, NOT ENOUGH VCF ENTRIES", pspos)
                    # re-initialize variables **from current SNP**
                    lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(
                        keys, l)
                    if test_leave(n, nmax) or not lgi:
                        break
                else:
                    n = check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto,
                                         sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb,
                                         "PASSED: ENOUGH VCF ENTRIES IN PREDEFINED WINDOW",
                                         "DISCARDED: NOT ENOUGH VCF ENTRIES IN PREDEFINED WINDOW", pspos)
                    lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(
                        keys, l)
                    if test_leave(n, nmax) or not lgi:
                        break
                    continue    # restart loop if predefined windows

            # 3.5) check SNP position < window start    (only for predefined
            # windows)
            elif int(spos) < int(wsta):
                if not bed:
                    sys.exit(
                        "BUG: vcf entry position < window start with sliding window!")
                l, pspos = update_SNP(f, spos)
                if not l:
                    break
                continue

            # 3.6) If allright, record allele info
            vgts = [ele[x] for x in inds]
            QUAL = float(ele[5])
            FORMAT = ele[8].split(":")
            # WARNING TODO MODIFY: HORRIBLE HACK, SPLIT THE ALT IF THERE IS MORE THAN ONE ALLELE
            # alleles = ele[3:5]
            alleles = [item for sublist in ele[3:5] for item in sublist.split(",")]
            seq, nbad = rec_allele(vgts, seq, keys, alleles, nbad, FORMAT, QUAL)
            # 3.7) check if enough data for all individual in the current window
            # i.e. no more than 'miss' missing SNPs
            if not bed:
                for j, nb in enumerate(nbad):
                    if nb > mmiss:
                        bad[j] = True
                if True in bad:
                    print_info("DISCARDED: TOO MUCH MISSING DATA PER INDIVIDUAL", wscaf, wsta, wsto, lgi,
                               sscaf, nbad, bad, pspos)
                    # thus change SNP and window!
                    l, pspos = update_SNP(f, spos)
                    if not l:
                        break
                    lgi, seq, nbad, bad, wscaf, wsta, wsto, wsize = initialize_win(
                        keys, l)
                    if not lgi:
                        break
                    continue

            # 3.8) increment number of analysed lines and proceed acccordingly.
            l, pspos = update_SNP(f, spos)
            if not l:
                break
            lgi += 1
        else:
            sys.exit("ERROR: problem vcf format line:\n%s" % l)
    return (lgi, minsnp, wscaf, wsta, wsto, lgi, sscaf, nbad, bad, n, keys,
            seq, outg, fnw, fout, vb, pspos, ele)


# SCRIPT
f = open_fil(vcf)
with gzip.open(nwk, 'wb') as fnw:
    with gzip.open(outc, 'wb') as fout:
        fout.write("scaffold\tstart\tend\tmid\tlength\tsites\n")
        l = f.readline()    # read vcf very first line
        # I) process vcf
        lgi, minsnp, wscaf, wsta, wsto, lgi, sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb, pspos, ele = main(
            l, indivs)

        # II) check last windows if enough snps analysed, then write the temp
        # fasta file
        if not bed:
            n = check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto,
                                 sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb,
                                 "LAST WINDOW PASSED: ENOUGH VCF ENTRIES IN SLIDING WINDOW",
                                 "LAST WINDOW DISCARDED: NOT ENOUGH VCF ENTRIES IN SLIDING WINDOW", pspos)
        elif lgi:
            n = check_lgi_minsnp(lgi, minsnp, wscaf, wsta, wsto,
                                 sscaf, nbad, bad, n, keys, seq, outg, fnw, fout, vb,
                                 "LAST WINDOW PASSED: ENOUGH VCF ENTRIES IN PREDEFINED WINDOW",
                                 "LAST WINDOW DISCARDED: NOT ENOUGH VCF ENTRIES IN PREDEFINED WINDOW", pspos)

f.close()
if bed:
    bed.close()
