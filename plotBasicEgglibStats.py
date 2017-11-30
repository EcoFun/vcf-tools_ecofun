#! /usr/bin/env python3

# ### WARNING : ONLY A COPY OF THE PLOTING PART OF computeBasicEgglibStats

import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import copy
import argparse

parser = argparse.ArgumentParser(description="""
Take the output of computeBasicEgglibStats and generate plots
""")
parser.add_argument('--statFiles', metavar='<statFiles>', nargs='*', help='The files containing all the statistics')
parser.add_argument('--outFolder', metavar='<outFolder>', type=str, help='The folder where (pdf) results will be written')
parser.add_argument('--subFile', metavar='<subFile>', type=str, help='A file containing a list of genes to be analyzed separatly')
args = parser.parse_args()
subFile = args.subFile
subList = [e.split(".")[0] for e in open(subFile).readlines()]
subList = list(set(subList))


def generateDictionaries(inFile, subList):
    dsub = dict()
    dother = dict()
    f = open(inFile)
    header = f.readline().strip().split("\t")
    for line in f.readlines():
        lineT = line.strip().split("\t")
        name = lineT[0].split(".")[0]
        for i in range(1, len(lineT)):
            val = lineT[i]
            stat = header[i]
            if name in subList:
                try:
                    dsub[stat][name] = cast(val)
                except:
                    dsub[stat] = dict()
                    dsub[stat][name] = cast(val)
            else:
                try:
                    dother[stat][name] = cast(val)
                except:
                    dother[stat] = dict()
                    dother[stat][name] = cast(val)
    return (dother, dsub)


def cast(value):
    ret = None
    try:
        ret = int(value)
    except:
        try:
            ret = float(value)
        except:
            if value == "None":
                return None
            ret = value.strip("[]").split(",")
    return ret


def plotStats(statDictionaryT, fileOut, color1=None, normed=False, statDictionaryTPartial=None, color2=None):
    statDictionary = copy.deepcopy(statDictionaryT)
    with PdfPages(fileOut) as pdf:
        for statName, statGenesD in sorted(statDictionary.items()):
            plt.figure()
            plt.title(statName)
            statGenes = cleanStats(statGenesD, 1e10)
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


def cleanStats(statsDic, limit):
    retList = list()
    for gene, value in statsDic.items():
        try:
            if value < -limit or value > limit:
                retList.append(None)
            else:
                retList.append(value)
        except:
            retList.append(None)
    return retList


def crap():
    if figOutName:
        outFile = os.path.join(outFolder, figOutName + ".pdf")
        if subList:
            plotStats(stats2[0], [stats2[1]], outFile, color1='yellow', statDictionaryTPartial=stats2[2], otherStatsTPartial=[stats2[3]], color2='red')
        else:
            plotStats(stats2[0], [stats2[1]], outFile, 'yellow')
    print(stats[-1], file=open(os.path.join(outFolder, figOutName + ".mk.res"), "w"))


for statFile in args.statFiles:
    a, b = generateDictionaries(statFile, subList)
    plotStats(a, "test" + statFile.split(".")[0] + ".pdf", color1='yellow', statDictionaryTPartial=b, color2='red')
