#!/usr/bin/env python

from __future__ import print_function
import glob
import egglib
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tpayenFunctions import addDic, encapsulate
import math
import copy
import numpy
import sys
import itertools

fastas = glob.glob("/home/tibo/09_gene_sequences/multipleAnalysis.filtered.genes/*fasta")


statisticsList = ['Aing', 'Aotg', 'As', 'Asd', 'Atot', 'B', 'Ch', 'ChE', 'D', 'Da', 'Deta', 'Dfl', 'Dj', 'Dstar', 'Dxy', 'E', 'F', 'Fs', 'Fstar', 'Gst', 'Gste', 'Hns', 'Hsd', 'Hst', 'K', 'Ke', 'Pi', 'Q', 'R', 'R2', 'R2E', 'R3', 'R3E', 'R4', 'R4E', 'Rintervals', 'Rmin', 'RminL', 'S', 'So', 'Ss', 'Sso', 'WCst', 'Z*nS', 'Z*nS*', 'ZZ', 'Za', 'ZnS', 'eta', 'etao', 'lseff', 'lseffo', 'nM', 'nPairs', 'nPairsAdj', 'ns_site', 'nseff', 'nseffo', 'nsingld', 'nsmax', 'nsmaxo', 'numFxA', 'numFxA*', 'numFxD', 'numFxD*', 'numShA', 'numShA*', 'numShP', 'numShP*', 'numSp', 'numSp*', 'numSpd', 'numSpd*', 'pM', 'rD', 'singl', 'singl_o', 'sites', 'sites_o', 'thetaH', 'thetaIAM', 'thetaL', 'thetaPi', 'thetaSMM', 'thetaW']

#statisticsList = ['Aing', 'Aotg', 'As', 'Asd', 'Atot', 'B', 'Ch', 'ChE', 'D', 'Da', 'Deta', 'Dfl', 'Dj', 'Dstar', 'Dxy', 'E', 'F', 'Fs', 'Fst', 'Fstar', 'Gst', 'Gste', 'Hns', 'Hsd', 'Hst', 'K', 'Ke', 'Kst', 'Pi', 'Q', 'R', 'R2', 'R2E', 'R3', 'R3E', 'R4', 'R4E', 'Rintervals', 'Rmin', 'RminL', 'S', 'Snn', 'So', 'Ss', 'Sso', 'WCst', 'Z*nS', 'Z*nS*', 'ZZ', 'Za', 'ZnS', 'eta', 'etao', 'lseff', 'lseffo', 'nM', 'nPairs', 'nPairsAdj', 'ns_site', 'nseff', 'nseffo', 'nsingld', 'nsmax', 'nsmaxo', 'numFxA', 'numFxA*', 'numFxD', 'numFxD*', 'numShA', 'numShA*', 'numShP', 'numShP*', 'numSp', 'numSp*', 'numSpd', 'numSpd*', 'pM', 'rD', 'singl', 'singl_o', 'sites', 'sites_o', 'thetaH', 'thetaIAM', 'thetaL', 'thetaPi', 'thetaSMM', 'thetaW']

namesFile = "/home/tibo/genesAnalysis/species.pop"

structureModification = {}  # {100: [300, 600, 700, 1000, 1100, 1200], 200: [500, 800, 900]}
toDel = [100, 200, 300, 500, 600, 700, 800, 900, 1000, 1100, 1200]  # [100, 200, 300, 500, 600, 700, 800, 900, 1000, 1100, 1200]
subList = [e.split(".")[0] for e in open("/home/tibo/08_gene_predictions/effectors/Vina_EUB04_domEU.v5pacbio.withExtrinsic.effectors.list").readlines() if e[0] != "#" and "Effector" in e]
subList = list(set(subList))

outFolder = "/home/tibo/genesAnalysis/egglibStats"


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


def computeStatsGenes(fastas, statsList, namesFile=None, modif=None, subList=None):
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
        nameGene = os.path.split(fa)[1].split(".")[0]
        if namesFile:
            try:
                aln = egglib.io.from_fasta(changeFastaFormat(fa, namesFile), string=True, groups=True, cls=egglib.Align)
                struct = egglib.stats.get_structure(aln, lvl_clust=0, lvl_pop=1, lvl_indiv=2)
                if modif:
                    struct = modifyStruct(struct, modif)
                # TODO: structD = struct.as_dict()
                # TODO: print("####")
                # TODO: print(structD[0])
                # TODO: print(len(structD[1]))
                tmp = cs.process_align(aln, struct=struct)
            except ValueError, e:
                print(fa, file=sys.stderr)
                count += 1
                print(count, file=sys.stderr)
                print(e, file=sys.stderr)
                continue
            except EgglibSizeException:
                    return None
        else:
            aln = egglib.io.from_fasta(fa)
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
        try:
            if codstat.num_codons_eff > 20:
                dnds = (codstat.num_pol_NS / codstat.num_sites_NS) / (codstat.num_pol_S / codstat.num_sites_S)
            else:
                dnds = None
        except ZeroDivisionError:
            print(codstat.num_codons_eff, file=sys.stderr)
            dnds = None
        dndsList.append(dnds)
        if subList and nameGene in subList:
            dndsSubList.append(dnds)
    if subList:
        return (stats, dndsList, statsSubList, dndsSubList)
    else:
        return (stats, dndsList)


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


def changeFastaFormat(fastaFile, namesFile):
    out = ""
    assocDic = dict()
    for nameLine in open(namesFile).readlines():
        nameTab = nameLine.strip().split()
        assocDic[nameTab[0]] = ">" + nameTab[0] + " @" + ",".join(nameTab[1:])
    for line in open(fastaFile).readlines():
        if line[0] != ">":
            out += line
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
        del dictionary[cluster][population]
    return dictionary


def main(fastas, statsList, namesFile=None, subList=None, structureModification={}, toDel=[], figOutName=None):
    if figOutName:
        print(figOutName)
    stats = computeStatsGenes(fastas, statsList, namesFile, (structureModification, toDel), subList)
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


if __name__ == "__main__":
    # name, toDelList = ("500vs1200", [100, 200, 300, 600, 700, 800, 900, 1000, 1100])
    # main(fastas, statisticsList, namesFile, subList=subList, structureModification=structureModification, toDel=toDelList, figOutName=name)
    for name, toDelList in createAllTuple(toDel, 1).items():
        print(name, file=sys.stderr)
        if os.path.exists(os.path.join(outFolder, name + ".pdf")):
            continue
        else:
            main(fastas, statisticsList, namesFile, subList=subList, structureModification=structureModification, toDel=toDelList, figOutName=name)
