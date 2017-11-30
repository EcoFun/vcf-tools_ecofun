#!/usr/bin/env python

from __future__ import print_function
import egglib
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tpayenFunctions import addDic, encapsulate, loadReferenceFasta
import copy
import numpy
import sys
import itertools
import argparse

parser = argparse.ArgumentParser(description="""
Create statistics for multi-fastas
""")
parser.add_argument('--fastas', metavar='<fastas>', nargs='*', help='The multi-fasta list of genes  that will be processed')
parser.add_argument('--namesFile', metavar='<namesFile>', type=str, help='The associative file between individuals and populations')
parser.add_argument('--outFolder', metavar='<outFolder>', type=str, help='The folder where (pdf) results will be written')
parser.add_argument('--subFile', metavar='<subFile>', type=str, help='A file containing a list of genes to be analyzed separatly')
args = parser.parse_args()

fastas = [os.path.abspath(e) for e in args.fastas]
namesFile = os.path.abspath(args.namesFile)
outFolder = os.path.abspath(args.outFolder)
subFile = os.path.abspath(args.subFile)

NB_POP_INCLUDED = 2


#### WARNING: WITH EGGLIB 3.0.0b12: Fst results in a segfault
statisticsList = ['Aing', 'Aotg', 'As', 'Asd', 'Atot', 'B', 'Ch', 'ChE', 'D', 'Da', 'Deta', 'Dfl', 'Dj', 'Dstar', 'Dxy', 'E', 'F', 'Fs', 'Fstar', 'Gst', 'Gste', 'Hns', 'Hsd', 'Hst', 'K', 'Ke', 'Pi', 'Q', 'R', 'R2', 'R2E', 'R3', 'R3E', 'R4', 'R4E', 'Rintervals', 'Rmin', 'RminL', 'S', 'So', 'Ss', 'Sso', 'WCst', 'Z*nS', 'Z*nS*', 'ZZ', 'Za', 'ZnS', 'eta', 'etao', 'lseff', 'lseffo', 'nM', 'nPairs', 'nPairsAdj', 'ns_site', 'nseff', 'nseffo', 'nsingld', 'nsmax', 'nsmaxo', 'numFxA', 'numFxA*', 'numFxD', 'numFxD*', 'numShA', 'numShA*', 'numShP', 'numShP*', 'numSp', 'numSp*', 'numSpd', 'numSpd*', 'pM', 'rD', 'singl', 'singl_o', 'sites', 'sites_o', 'thetaH', 'thetaIAM', 'thetaL', 'thetaPi', 'thetaSMM', 'thetaW']

structureModification = {}  # EXAMPLE {100: [300, 600, 700, 1000, 1100, 1200], 200: [500, 800, 900]}
subList = [e.split(".")[0] for e in open(subFile).readlines()]
subList = list(set(subList))
individualsToDel = ["Vina-2269-pyrac-trim", "Vina-2504-loq-trim", "Vina-2263-loq-trim", "Vina-2474-CAP-trim", "Vina-2225-CAM-trim", "Vina-2226-CAM-trim", "Vina-2458-CAM-trim", "Vina-2478-syl-trim"]


def computeMK(align, struct):
    # Computed if NB_POP_INCLUDED = 1
    # Ds: the number of synonymous substitutions per gene
    # Dn: the number of non-synonymous substitutions per gene
    # Ps: the number of synonymous polymorphisms per gene
    # Pn: the number of non-synonymous polymorphisms per gene
    # alpha: 1 - (Ds * Pn) / (Dn * Ps)
    Ds = 0
    Dn = 0
    Ps = 0
    Pn = 0
    cs = egglib.stats.ComputeStats()
    cs.add_stats('numFxD')
    cs.add_stats('S')
    cdiv = egglib.stats.CodingDiversity(align)
    cdiv_S = cdiv.mk_align_S()
    cdiv_NS = cdiv.mk_align_NS()
    resS = cs.process_align(cdiv_S, struct=struct, filtr=egglib.stats.filter_codon)
    resNS = cs.process_align(cdiv_NS, struct=struct, filtr=egglib.stats.filter_codon)
    Ds = resS['numFxD']
    Dn = resNS['numFxD']
    Ps = resS['S']
    Pn = resNS['S']
    try:
        alpha = 1 - (Ds * Pn) / (float(Dn) * Ps)
    except ZeroDivisionError:
        alpha = None
    except TypeError:
        alpha = None
    return {"Ds": Ds, "Dn": Dn, "Ps": Ps, "Pn": Pn, "alpha": alpha}


def computeStatsGenes(fastas, statsList, namesFile=None, modif=None, subList=None, individualsToDelete=None, deletedIndividualsStructure=None):
    dicoPNPSDNDS = dict()
    cs = egglib.stats.ComputeStats()
    for s in statsList:
        cs.add_stats(s)
    stats = dict()
    dndsList = list()
    statsSubList = dict()
    dndsSubList = list()
    errors = 0
    count = 0
    for fa in fastas:
        if individualsToDelete:
            fastaString = deleteIndividuals(fa, individualsToDelete, deletedIndividualsStructure)
        else:
            fastaString = open(fa).read()
        nameGene = os.path.split(fa)[1].split(".")[0]
        if namesFile:
            try:
                print("In 1", file=sys.stderr)
                aln = egglib.io.from_fasta(changeFastaFormat(fastaString, namesFile), string=True, groups=True, cls=egglib.Align)
                print("In 2", file=sys.stderr)
                struct = egglib.stats.get_structure(aln, lvl_clust=0, lvl_pop=1, lvl_indiv=2)
                print("In 3", file=sys.stderr)
                if modif:
                    struct = modifyStruct(struct, modif)
                print("In 4", file=sys.stderr)
                tmp = cs.process_align(aln, struct=struct)
                print("In 5", file=sys.stderr)
            except ValueError, e:
                print(fa, file=sys.stderr)
                count += 1
                print(count, file=sys.stderr)
                print(e, file=sys.stderr)
                tmp = dict()
                for e in statsList:
                    tmp[e] = None
            except EgglibSizeException:
                tmp = dict()
                for e in statsList:
                    tmp[e] = None
            try:
                MK = computeMK(aln, struct)
            except:
                MK = {}
            dicoPNPSDNDS[nameGene] = MK
        else:
            aln = egglib.io.from_fasta(fastaString, string=True)
            tmp = cs.process_align(aln)
        try:
            stats = addDic(stats, encapsulate(tmp))
        except:
            stats = encapsulate(tmp.copy())
        if subList and nameGene in subList:
            try:
                statsSubList = addDic(statsSubList, encapsulate(tmp))
            except:
                statsSubList = tmp.copy()
        codstat = egglib.stats.CodingDiversity()
        try:
            codstat.process(aln)
        except TypeError:
            errors += 1
        except Exception, e:
            errors += 1
            print(e, file=sys.stderr)
        try:
            if codstat.num_codons_eff > 20:
                dnds = (codstat.num_pol_NS / codstat.num_sites_NS) / (codstat.num_pol_S / codstat.num_sites_S)
            else:
                dnds = None
        except ZeroDivisionError:
            print(codstat.num_codons_eff, file=sys.stderr)
            dnds = None
        except e:
            print(e, file=sys.stderr)
            dnds = None
        dndsList.append(dnds)
        if subList and nameGene in subList:
            dndsSubList.append(dnds)
    if subList:
        return (stats, dndsList, statsSubList, dndsSubList, dicoPNPSDNDS)
    else:
        return (stats, dndsList, dicoPNPSDNDS)


def plotStats(statDictionaryT, otherStatsT, fileOut, color1=None, normed=False, statDictionaryTPartial=None, otherStatsTPartial=None, color2=None):
    statDictionary = copy.deepcopy(statDictionaryT)
    otherStats = copy.deepcopy(otherStatsT)
    with PdfPages(fileOut) as pdf:
        for statName, statGenes in statDictionary.items():
            plt.figure()
            plt.title(statName)
            statGenes = cleanStats(statGenes, 1e10)
            try:
                data = [round(e, 2) for e in statGenes if e]
                try:
                    min(data)
                    plt.hist(data, log=False, color=color1, normed=normed, alpha=0.5)
                except ValueError:
                    pass
                if statDictionaryTPartial:
                    statGenesPartial = cleanStats(statDictionaryTPartial[statName], 1e10)
                    data = [round(e, 2) for e in statGenesPartial if e]
                    try:
                        min(data)
                        plt.hist(data, log=False, color=color2, normed=normed, alpha=0.5)
                    except ValueError:
                        pass
            except TypeError:
                pass
            pdf.savefig()
            plt.close()
        for stat in otherStats:
            plt.hist([e for e in stat if e], bins=100, log=False, color=color1, normed=normed)
            if otherStatsTPartial:
                statPartial = otherStatsTPartial.pop(0)
                plt.hist([e for e in statPartial if e], bins=100, log=False, color=color2, normed=normed)
            plt.title('dN/dS')
            pdf.savefig()
            plt.close()


def cleanStats(listOfStats, limit):
    for stat in range(len(listOfStats)):
        if listOfStats[stat] < -limit or listOfStats[stat] > limit:
            listOfStats[stat] = None
    return listOfStats


def writeStats(statDictionaryT, otherStatsT, fastas):
    statDictionary = copy.deepcopy(statDictionaryT)
    otherStats = copy.deepcopy(otherStatsT)
    print("\t" + "\t".join(statDictionary.keys()))
    means = list()
    standardDeviation = list()
    for v in statDictionary.values():
        try:
            means.append(str(numpy.mean([e for e in v if e])))
        except:
            means.append("None")
        try:
            standardDeviation.append(str(numpy.std([e for e in v if e])))
        except:
            standardDeviation.append("None")
    print("Means\t" + "\t".join(means))
    print("Std\t" + "\t".join(standardDeviation))
    for fa in fastas:
        line = os.path.split(fa)[1] + "\t"
        for val in statDictionary.values():
            line += str(val.pop(0)) + "\t"
        for stat in otherStats:
            line += str(stat.pop(0)) + "\t"
        print(line)


def createAllTuple(listOfPopulations, sizeOfTuple):
    outs = dict()
    tupleKept = set([tuple(sorted(e)) for e in itertools.permutations(listOfPopulations, sizeOfTuple)])
    for kept in tupleKept:
        listToRemove = copy.deepcopy(listOfPopulations)
        name = "vs".join(tuple([str(e) for e in kept]))
        for elem in kept:
            listToRemove.remove(elem)
        outs[name] = listToRemove
    return outs


def modifyStruct(struct, toChange):
    modif = toChange[0]
    populationList = toChange[1]
    dicIngroup, dicOutgroup = struct.as_dict()
    deletePopulation(dicIngroup, 1000, populationList)
    for k, v in modif.items():
        for elem in v:
            moveGroup(dicIngroup, 1, elem, k)
    for cluster in dicIngroup.values():
        for population in cluster.values():
            if len(population) < 3:
                raise EgglibSizeException
    struct.make_structure(dicIngroup, dicOutgroup)
    return struct


class EgglibSizeException(Exception):
    """All population must be larger than one"""


def deleteIndividuals(fastaFile, individualsToDelete, secondaryList=None):
    fastaDic = loadReferenceFasta(fastaFile)
    out = ""
    numberDeleted = 0
    for k, v in fastaDic.items():
        if k in individualsToDelete or (secondaryList and k not in secondaryList):
            numberDeleted += 1
        else:
            try:
                out += ">" + k + "\n" + v + "\n"
            except:
                print("Check if " + fastaFile + " is missing", file=sys.stderr)
    if numberDeleted != len(individualsToDelete):
        print("In: ", fastaFile, " there were ", numberDeleted, " individuals deleted and ", len(individualsToDelete), " were expected", file=sys.stderr)
    return out


def changeFastaFormat(fastaString, namesFile):
    out = ""
    assocDic = dict()
    for nameLine in open(namesFile).readlines():
        nameTab = nameLine.strip().split()
        assocDic[nameTab[0]] = ">" + nameTab[0] + " @" + ",".join(nameTab[1:])
    for line in fastaString.split("\n"):
        if line == "":
            continue
        if line[0] != ">":
            out += line + "\n"
        else:
            name = line.strip(">").strip()
            out += assocDic[name] + "\n"
    return out


# Always 3 lvl
def moveGroup(dictionary, lvl, oldValue, newValue):
    if lvl > 2:
        raise Exception("Expect the level to be in [0-2]")
    if lvl == 0:
        try:
            tempD = dictionary.pop(oldValue)
        except:
            raise Exception("The old value must exist")
        if newValue:
            dictionary[newValue].update(tempD)
    if lvl == 1:
        success = False
        for k in dictionary.keys():
            try:
                tempD = dictionary[k].pop(oldValue)
                success = True
            except:
                continue
            if newValue:
                dictionary[k][newValue].update(tempD)
        if not success:
            raise Exception("The old value must exist")
    if lvl == 3:
        success = False
        for k1 in dictionary.keys():
            for k2 in dictionary[k1].keys():
                try:
                    tempD = dictionary[k1][k2].pop(oldValue)
                    success = True
                except:
                    continue
                if newValue:
                    dictionary[k1][k2][newValue].update(tempD)
        if not success:
            raise Exception("The old value must exist")
    return dictionary


def movePopulationCluster(dictionary, oldCluster, newCluster, populationList):
    for population in populationList:
        if newCluster not in dictionary.keys():
            dictionary[newCluster] = dict()
        dictionary[newCluster][population] = dictionary[oldCluster].pop(population)
    return dictionary


def deletePopulation(dictionary, cluster, populationList):
    for population in populationList:
        try:
            del dictionary[cluster][population]
        except:
            print(str(cluster) + "." + str(population) + " is already non existant", file=sys.stderr)
    return dictionary


def main(fastas, statsList, namesFile=None, subList=None, structureModification={}, toDel=[], figOutName=None, individualsToDelete=None, indivDeletedStructure=None):
    if figOutName:
        print(figOutName)
    stats = computeStatsGenes(fastas, statsList, namesFile, (structureModification, toDel), subList=subList, individualsToDelete=individualsToDelete, deletedIndividualsStructure=indivDeletedStructure)
    if not stats:
        print("Stats are missing", file=sys.stderr)
        return None
    stats2 = stats
    writeStats(stats[0], [stats[1]], fastas)
    if figOutName:
        outFile = os.path.join(outFolder, figOutName + ".pdf")
        if subList:
            plotStats(stats2[0], [stats2[1]], outFile, color1='yellow', statDictionaryTPartial=stats2[2], otherStatsTPartial=[stats2[3]], color2='red')
        else:
            plotStats(stats2[0], [stats2[1]], outFile, 'yellow')
    print(stats[-1], file=open(os.path.join(outFolder, figOutName + ".mk.res"), "w"))

    
if __name__ == "__main__":
    toDel = sorted(list(set([int(e.split()[2]) for e in open(namesFile).readlines() if len(e.split()) > 2])))
    for name, toDelList in createAllTuple(toDel, NB_POP_INCLUDED).items():
        print(name, file=sys.stderr)
        if os.path.exists(os.path.join(outFolder, name + ".pdf")):
            continue
        else:
            main(fastas, statisticsList, namesFile, subList=subList, structureModification=structureModification, toDel=toDelList, figOutName=name, individualsToDelete=individualsToDel, indivDeletedStructure=[e.split()[0] for e in open(namesFile).readlines()])

