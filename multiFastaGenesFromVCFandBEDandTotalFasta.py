#! /usr/bin/env python

from __future__ import print_function
import argparse
import vcf
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description='Create one multifasta by gene in the BED file')

parser.add_argument('input_vcf', help='vcf to process. BEWARE: the vcf must be sorted per coordinates!', nargs=1)
parser.add_argument('input_annotation', help='annotation to process (in gff3 for now). BEWARE: the annotation must be sorted per coordinates!', nargs=1)
parser.add_argument('input_fasta', help='Fasta of the complete reference genome.', nargs=1)
# TODO: For now not possible to change
parser.add_argument('output_fasta_folder', help='output folder. If not given write all the fasta in place (For now not possible to change)', nargs=1, default=".")
args = parser.parse_args()
nameReferenceGenome = "Vina-EUB04-domEU"


def getAnnotation(annotLine):
    ret = getAnnotation_gff(annotLine)
    return ret


def getAnnotation_gff(annotationLine):
    if "CDS" not in annotationLine:
        return False
    annotationInterTable = annotationLine.strip().split("\t")
    if len(annotationInterTable) < 5:
        return False
    annotationInterSca = annotationInterTable[0]
    annotationInterStart = int(annotationInterTable[3])
    annotationInterStop = int(annotationInterTable[4])
    annotationGene = annotationInterTable[8].split(";")[1].split()[1].strip("\"")  # THIS PART IS PROBABLY HIGHLY SPECIFIC
    annotationSens = annotationInterTable[6]   # Should be the same for all CDS in the same gene
    return (annotationInterSca, annotationInterStart, annotationInterStop, annotationGene, annotationSens)


def writeAbsentVCF(sampleNames, length=1):
    tempDic = dict()
    for name in sampleNames:
        tempDic[name] = "N" * length
    return tempDic


def addDic(dic1, dic2):
    tmpDic = dict()
    if (len(dic1) != (len(dic2) + 1)):
        raise Exception("Dcitionnaries ar not of the same size")
    else:
        for k, v in dic1.items():
            try:
                tmpDic[k] = v + dic2[k]
            except KeyError:
                tmpDic[k] = v
            except:
                raise Exception
    return tmpDic


def reverse_complement(dic):
    tmpDic = dict()
    for k, v in dic.items():
        tmpDic[k] = str(Seq(v).reverse_complement())
    return tmpDic


def increaseVCF(vcf_reader, annotSca, annotStart):
    try:
        vcfInfo = vcf_reader.next()
        while (vcfInfo.CHROM <= annotSca and vcfInfo.POS < annotStart):
            vcfInfo = vcf_reader.next()
        return vcfInfo
    except:
        raise Exception("VCF Finished")


def generatePosition(mutation):
    tmpDic = dict()
    for s in mutation.samples:
        try:
            if len(s.gt_bases) > 1:
                print(s.gt_bases)
            allele = s.gt_bases
        except:
            allele = "N"
        tmpDic[s.sample] = allele  # TODO CHECK QUALITY
    return tmpDic


def initGeneSequence(samples):
    tmpDic = dict()
    for sample in samples:
        tmpDic[sample] = ""
    return tmpDic


# PROCESS
def main(vcf_reader, annotationHandle, refDic):
    sampleNames = vcf_reader.samples
    # Annotation
    presentGeneName = ""
    presentgeneSens = ""
    presentGeneSequence = initGeneSequence(sampleNames)
    # Mutation
    mutationPosition = vcf_reader.next()
    # for position in vcf_reader
    mutationScaffold = mutationPosition.CHROM
    mutationPositionNumber = mutationPosition.POS
    #
    # We read all the annotation file
    #
    for annotationInterval in annotationHandle:
        # We only want the CDS intervals
        annotInfo = getAnnotation(annotationInterval)
        if annotInfo:
            annotationIntervalSca = annotInfo[0]
            annotationIntervalStart = annotInfo[1]
            annotationIntervalStop = annotInfo[2]
            annotationGene = annotInfo[3]
            annotationSens = annotInfo[4]
        # If annotInfo is False we go to the next information (not a CDS)
        else:
            continue
        # We write the previous gene
        if (annotationGene != presentGeneName):
            if presentGeneName != "":
                if presentgeneSens == "-":
                    presentGeneSequence = reverse_complement(presentGeneSequence)
                print(presentGeneName + " is generated")
                # TODO Permit to change the output folder
                print("\n".join([">" + k + "\n" + v + "\n" for k, v in sorted(presentGeneSequence.items(), key=lambda e: e[0])]), file=open(presentGeneName + ".multi.fasta", "w"))
                # We put the new informations in place
                presentGeneName = annotationGene
                presentgeneSens = annotationSens
                presentGeneSequence = initGeneSequence(sampleNames)
            else:
                presentGeneName = annotationGene
        refSeq = refDic[annotationIntervalSca][annotationIntervalStart - 1:annotationIntervalStop]
        try:
            presentGeneSequence[nameReferenceGenome] += refSeq
        except:
            presentGeneSequence[nameReferenceGenome] = refSeq
        # If the annotationInterSca is after the scaffold of the vcf go to the next position of the vcf
        if (annotationIntervalSca > mutationScaffold):
            try:
                mutationPosition = increaseVCF(vcf_reader, annotationIntervalSca, annotationIntervalStart)
            except:
                return 0
            mutationScaffold = mutationPosition.CHROM
            mutationPositionNumber = mutationPosition.POS
        # We have gone to the next interesting vcf position

        # If the bedInterSca is before the scaffold of the vcf go to the next position of the bed
        if (annotationIntervalSca < mutationScaffold):
            # ALL THE ANNOTATION HAVE NO CORRESPONDANT IN VCF, ENTIRE ANNOTATION WITH ONLY N FOR ALL BUT REFERENCE
            print("Add N scaffold")  # TODO REMOVE
            presentGeneSequence = addDic(presentGeneSequence, writeAbsentVCF(sampleNames, annotationIntervalStop - annotationIntervalStart + 1))
            continue

        if (annotationIntervalStart > mutationPositionNumber):
            try:
                mutationPosition = increaseVCF(vcf_reader, annotationIntervalSca, annotationIntervalStart)
            except:
                return 0
            mutationScaffold = mutationPosition.CHROM
            mutationPositionNumber = mutationPosition.POS
        # We have gone to the next interesting vcf position

        # The annotation interval is completely before the vcf position, we go to next group of position in annotation
        if (annotationIntervalStop < mutationPositionNumber):
            # ALL THE ANNOTATION HAVE NO CORRESPONDANT IN VCF, ENTIRE ANNOTATION WITH ONLY N FOR ALL BUT REFERENCE
            print("Add N batch")  # TODO REMOVE
            presentGeneSequence = addDic(presentGeneSequence, writeAbsentVCF(sampleNames, annotationIntervalStop - annotationIntervalStart + 1))
            continue
        # The annotation interval is after the vcf position, we go to next position of vcf

        # The vcf position is in the BED interval
        if (mutationPositionNumber > annotationIntervalStart):
            presentGeneSequence = addDic(presentGeneSequence, writeAbsentVCF(sampleNames, mutationPositionNumber - annotationIntervalStart))

        positionPointer = mutationPositionNumber
        while (mutationScaffold == annotationIntervalSca and mutationPositionNumber <= annotationIntervalStop):
            # Check if the mutations follow one another, if it do, generate and go to the next
            if mutationPositionNumber == positionPointer:
                d = generatePosition(mutationPosition)
                presentGeneSequence = addDic(presentGeneSequence, d)
                try:
                    mutationPosition = increaseVCF(vcf_reader, annotationIntervalSca, annotationIntervalStart)
                except:
                    return 0
                mutationScaffold = mutationPosition.CHROM
                mutationPositionNumber = mutationPosition.POS
                positionPointer += 1
            # If mutations don't follow directy the previous one add N
            else:
                presentGeneSequence = addDic(presentGeneSequence, writeAbsentVCF(sampleNames))
                positionPointer += 1
        while len(presentGeneSequence[nameReferenceGenome]) > len(presentGeneSequence.values()[0]):
            presentGeneSequence = addDic(presentGeneSequence, writeAbsentVCF(sampleNames))


def loadReferenceFasta(fastaFileName):
    tmpDic = dict()
    f = open(fastaFileName)
    chro = None
    seq = ""
    for l in f.readlines():
        if l[0] == ">":
            if chro:
                tmpDic[chro] = seq
                seq = ""
            chro = l.split()[0].strip(">")
        else:
            seq += l.strip()
    tmpDic[chro] = seq
    return tmpDic


# BEGIN: OPENING
vcfHandle = open(args.input_vcf[0])
vcf_reader = vcf.Reader(vcfHandle)
annotationHandle = open(args.input_annotation[0])
dicReference = loadReferenceFasta(args.input_fasta[0])

main(vcf_reader, annotationHandle, dicReference)

# END: CLOSING
vcfHandle.close()
annotationHandle.close()


# for position in vcf_reader:
#    scaffold = position.CHROM
#    positionNumber = position.POS
