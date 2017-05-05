#!/bin/env python2

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

import sys, argparse as ap, os, gzip, re

# 0.1) get options from commands lines
parser = ap.ArgumentParser(description='Filter snpeff annotated vcf to detect interesting genes')

parser.add_argument('group_list', help='text file of individuals included in the focal group (one individual per line).', nargs=1)
parser.add_argument('input_vcf', help='snpeff vcf to process.', nargs=1)
parser.add_argument('output_csv', help='Csv file detailling all features with high effect variants', nargs=1)
parser.add_argument('output_list', help='List of genes derived from the list of features with high effect variants', nargs=1)

parser.add_argument('-g', '--gzip', help='Gzip result file [false].', action='store_true')
parser.add_argument('-m', '--moderate', help='Include Moderate in alleles of strong effects [False].', action='store_true')

# 0.2) set options
def make_gz(fil):
	if fil[-2:] != "gz":
		fil = fil + ".gz"
		return(fil)
args = parser.parse_args()
fgrp = args.group_list[0]
vcf = args.input_vcf[0]
do_gzip = args.gzip
foutp = make_gz(args.output_csv[0]) if do_gzip else args.output_csv[0]
foutp2 = args.output_list[0]
imod = args.moderate

#### FUNCTIONS
def open_fil(fil):
	if fil[-2:]=="gz":
		f = gzip.open(fil)
	else:
		f = open(fil)
	return(f)

def getlist(fil):
	lis = []
	f = open_fil(fil)
	for l in f:
		lis.append(l.strip())
	return(lis)

def get_real_allele(bin_alleles, alleles):
	return([ alleles[int(x)] for x in bin_alleles ] )

def get_og_all(i_ingp, ll, alleles):
	ogpallnb = [ ll[x].split(":")[0] for x in iogp ]
	ogpall = get_real_allele(ogpallnb, alleles)
	return(ogpall)

def get_effect_allele(pattern, lis, sep="|", field=0, get_ann_patt = False):
	ann_patt = filter(lambda x: pattern in x, lis)
	eff_allele = [ x.split(sep)[field] for x in ann_patt ]
	if get_ann_patt:
		res = eff_allele, ann_patt
	else:
		res = eff_allele
	return(res)

def get_LOF_effect_allele(gn, pattern, lis, sep="|", field=0):
	ann_patt = filter(lambda x: gn in x, lis)
	ann_patt = filter(lambda x: pattern in x, lis)
	eff_allele = [ x.split(sep)[field] for x in ann_patt ]
	res = eff_allele
	return(res)

def get_stg_effect_alleles(info, imod):
		# HIGH alleles
	linfo = info.split(";")
	ann = filter(lambda x: "ANN=" in x, linfo)[0]
	lann = re.split("=|,", ann)[1:]
	ahigh, ann_high = get_effect_allele("|HIGH|", lann, get_ann_patt = True)
		# LOF alleles
	lof = filter(lambda x: "LOF=" in x, linfo)
	if lof:
		ann_lof = lof[0]
		glof = re.split("=|\|", lof[0])[2]
		alofs = get_LOF_effect_allele(glof, "|HIGH|", lann)
	else:
		alofs, ann_lof, glof = [], '', ''
	lstrg = set(ahigh + alofs)
	res = lstrg, ann_high, ann_lof, glof, ann
	return(res)

def get_method(ll):
	if ll[8]=='GT:GQ:DP:AD:RO:QR:AO:QA:GL':
		test_method = "freebayes"
	elif ll[8]== 'GT:AD:DP:GQ:PL':
		test_method = "GATK-UG"
	else:
		test_method = None
	return(test_method)

def initialize_fe_dict(fe, grp, candidates):
	candidates.update ({ 
		fe: {
			'test_indiv_multi': dict((i, False) for i in grp), 'GATK-UG': 0, 
			'test_fixed': False, 'test_multiple': False , 'test_LOF': False, 'freebayes': 0,
			'fixed': {}, 'multiple': {}, 'fixed_ct':0, 'multiple_ct': 0,
			'partial_multiple_ct': 0, 'GATK-UG_fixed': 0, 'freebayes_fixed': 0 }})
	return(candidates)

def update_res_dict(mydic, method, key, temp_dic):
	mydic[method] += 1
	if mydic[key] == {}:
		mydic[key] = { 
			'feature': tuple([fe]) , 'gene': fe.split(".")[0], 'POS': "", 'REF': "", 'ALT': "",
			'Detection_method': "", 'HIGHeffectAlleles': "", 'LOFgene': "", 'LOFfield': "",
			'AnnHIGHeffects': "", 'fullAnn': "", 'Individuals': ""
			}
	for k in temp_dic.keys():
		if mydic[key][k] == "":
			mydic[key][k] += temp_dic[k]
		else:
			mydic[key][k] += ";" + temp_dic[k]
	return(mydic)

def update_tests(key, mydic, grp, test_stg_gts_ingp, method, tlof):
	if key == 'fixed':
		mydic['test_fixed'] = True
		mydic['fixed_ct'] += 1
		if method == "freebayes":
			mydic['freebayes_fixed'] += 1
		elif method == "GATK-UG":
			mydic['GATK-UG_fixed'] += 1
	elif key == 'multiple' and mydic['test_multiple'] == False:
		tindiv = dict((i, j) for i,j in zip(grp, test_stg_gts_ingp))
		mydic['partial_multiple_ct'] += 1
		for k in mydic['test_indiv_multi'].keys():
			if not mydic['test_indiv_multi'][k] and tindiv[k]:
				mydic['test_indiv_multi'][k] = True
		if False not in mydic['test_indiv_multi'].values():
			mydic['test_multiple'] = True
			mydic['multiple_ct'] += 1
	if not mydic['test_LOF'] and tlof:
		mydic['test_LOF'] = True
	return(mydic)

#### SCRIPT
f = open_fil(vcf)
grp = getlist(fgrp)
candidates = {}
for l in f:
	l = l.strip()
	if l[0] == "#":
		#~outp.write(l + "\n")
		if l[0:2] == "#C":
			ll = l.split("\t")
			indiv = ll[9:]
			try:
				i_ingp = [ indiv.index(x) for x in grp ]
				i_outgp = filter(lambda x: x not in i_ingp, range(0,len(indiv)))
			except ValueError:
				sys.exit(
				"""One or several individuals of the 'group_list' are not in the vcf header. 
				Please check the spelling.""")
	else:
		# 1) keep SNPs with strong effects only (else continue) - no missing data allowed
			# 1.1) test the presence of LOF, NMD or HIGH effects
		ll = l.split("\t")
		info = ll[7]
					# test LOF
		tlof = "LOF=" in info
		tnmd = "NMD=" in info
		thigh = "|HIGH|" in info
		# implement imod just in case for the moment
		if imod:
			tmod = "|MODERATE|" in info
			stg_effect = [ tlof, tnmd, thigh, tmod ]	# presence of strong effect?
		else:
			stg_effect = [ tlof, tnmd, thigh ]
		
		gtsbin = [ x.split(":")[0] for x in ll[9:] ]
		missdata = "." in gtsbin
			# 1.2) test: strong effect and no missing data
		if True not in stg_effect or missdata:
			continue
		else:
		# 2) test that none of the outgroup have a strong effect allele (else continue)
			# 2.1) get alleles of strong effect and LOF effect
			stg_alleles, ann_high, ann_lof, gn_lof, ann = get_stg_effect_alleles(info, imod)

			# 2.2) test that no individual of the outgroup have strong effect SNPs
			alleles = [ll[3]] + ll[4].split(",")
			gts = get_real_allele(gtsbin, alleles)
			gts_outgp = [ gts[x] for x in i_outgp ]
			stg_in_outgp = True in [ x in stg_alleles for x in gts_outgp ]
			if stg_in_outgp:
				continue
			else:
				# 3) are the strong effect alleles present in all individuals of the focal group?
				gts_ingp = [ gts[x] for x in i_ingp ]
				test_stg_gts_ingp = [ x in stg_alleles for x in gts_ingp ]
				if True not in test_stg_gts_ingp:
					continue
				all_ingp_stg = False not in test_stg_gts_ingp
				features = set([ x.split("|")[6] for x in ann_high ])
				detection_method = get_method(ll)
				# print in a output file
				if not detection_method:
					print "WARNING: Detection method is not GATK nor freebayes"
				if all_ingp_stg:
					key = 'fixed'
				else: 
					key = 'multiple'
				# update result dict
				for fe in list(features):
					gd_indiv = 'all' if all_ingp_stg else ','.join([ x for x,y in zip(grp, test_stg_gts_ingp) if y is True ])
					temp_dic = { 
						'POS': ll[0]+":"+ll[1], 'REF': ll[3], 'ALT': ll[4], 
						'Detection_method': detection_method, 'LOFgene': gn_lof, 
						'HIGHeffectAlleles': ",".join(stg_alleles), 'LOFfield': ann_lof,
						'AnnHIGHeffects': ",".join(ann_high), 'fullAnn': ann,
						'Individuals': gd_indiv }
					# if feature not in the final dict, create fe dict
					if fe not in candidates.keys():
						candidates = initialize_fe_dict(fe, grp, candidates)
					# update
					candidates[fe] = update_res_dict(candidates[fe], detection_method, key, temp_dic)
					candidates[fe] = update_tests(key, candidates[fe], grp, test_stg_gts_ingp, detection_method, tlof)

f.close()

## final list of results
final = {}
for k in candidates.keys():
	if candidates[k]['test_fixed'] or candidates[k]['test_multiple']:
		final[k] = candidates[k]


Ncandi = len(candidates)
NgNcandi = len(set([x.split(".")[0] for x in candidates.keys()]))
Nfinal = len(final)
NgNfinal = len(set([x.split(".")[0] for x in final.keys()]))

SNPfix_stg = 0 ; SNPpol_stg = 0
NgNmultiSNP = 0
for k in final.keys():
	SNPfix_stg += final[k]['fixed_ct']
	if final[k]['fixed_ct']>1 and final[k]['test_fixed'] == True:
		NgNmultiSNP += 1
	SNPpol_stg += final[k]['partial_multiple_ct']


### print results
	# write variant effects
if do_gzip:
	outp = gzip.open(foutp, 'wb')
else: 
	outp = open(foutp, 'wb')
gns = [ int(re.split("g|\.", x)[1]) for x in final.keys()]
isgn = [i[0] for i in sorted(enumerate(gns), key=lambda x:x[1])]
keys = [ final.keys()[x] for x in isgn ]
header = "#Feature\tGene\tfixed_ct\tPOS\tDetection_method\tREF\tHIGHeffectAlleles\tLOFfield\tAnnHIGHeffects\tfullAnn"
outp.write(header + "\n")
dfgen = {}
for k in keys:
	dk = final[k]
	gn = k.split(".")[0]
	try:
		if not dfgen[gn]['LOF'] and dk['test_LOF']:
			dfgen[gn]['LOF'] == True 
		dfgen[gn]['ct_variants'] += dk['fixed_ct'] + dk['multiple_ct']
		dfgen[gn]['N_features'] += 1
		if dk['freebayes'] and not dfgen[gn]['freebayes']:
			dfgen[gn]['freebayes'] = 1
		if dk['GATK-UG'] and not dfgen[gn]['GATK-UG']:
			dfgen[gn]['GATK-UG'] = 1
	except KeyError:
		dfgen[gn] = { 'LOF': dk['test_LOF'], 'ct_variants': dk['fixed_ct'] + dk['multiple_ct'], 'N_features': 1 , 'freebayes': dk['freebayes'],
		'GATK-UG': dk['GATK-UG'] }
	lres = [ k, gn, str(dk['fixed_ct']), dk['fixed']['Detection_method'], 
				dk['fixed']['POS'], dk['fixed']['REF'],
				dk['fixed']['HIGHeffectAlleles'], dk['fixed']['LOFfield'], 
				dk['fixed']['AnnHIGHeffects'], dk['fixed']['fullAnn'] ]
	res = "\t".join(lres)
	outp.write(res + "\n")
outp.close()

	# write final list of genes
		# include if LOF effect, variant number count
outp2 = open(foutp2, 'wb')
gns = dfgen.keys()
isgn = [i[0] for i in sorted(enumerate(gns), key=lambda x:x[1])]
keys = [ gns[x] for x in isgn ]
outp2.write("gene\tLOF\tN_features\tN_variants\tDetection_methods" + "\n")
for k in keys:
	det_meth = ""
	if dfgen[k]['GATK-UG']:
		det_meth += 'GATK-UG'
		if dfgen[k]['freebayes']:
			det_meth += ',freebayes'
	elif dfgen[k]['freebayes']:
		det_meth += 'freebayes'
	else:
		print "problem"
	outp2.write("%s\t%s\t%s\t%s\t%s\n" % (k, dfgen[k]['LOF'], dfgen[k]['N_features'], dfgen[k]['ct_variants'], det_meth))
outp2.close()






