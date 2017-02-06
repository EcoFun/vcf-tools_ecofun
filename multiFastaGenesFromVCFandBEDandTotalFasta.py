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


def getAnnotation(annotLine):
    ret = getAnnotation_gff(annotLine)
    return ret


def getAnnotation_gff(annotationLine):
    if "CDS" not in annotationLine:
        return False
    annotationInterTable = annotationLine.strip().split("\t")
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
    if (sorted(dic1.keys()) != sorted(dic2.keys())):
        raise Exception
    else:
        for k, v in dic1.items():
            tmpDic[k] = v + dic2[k]
    return tmpDic


def reverse_complement(dic):
    tmpDic = dict()
    for k, v in dic.items():
        tmpDic[k] = str(Seq(v).reverse_complement())
    return tmpDic


def increaseVCF(vcf_reader, annotSca, annotStart):
    vcfInfo = vcf_reader.next()
    while (vcfInfo.CHROM <= annotSca and vcfInfo.POS < annotStart):
        vcfInfo = vcf_reader.next()
    return vcfInfo


def generatePosition(mutation):
    tmpDic = dict()
    for s in mutation.samples:
        try:
            allele = mutation[s.e.gt_alleles[0]]
        except:
            allele = "N"
        tmpDic[s.sample] = allele  # TODO CHECK QUALITY
    return tmpDic


# PROCESS
def main(vcf_reader, annotationHandle):
    sampleNames = vcf_reader.samples
    # Annotation
    presentGeneName = ""
    presentgeneSens = ""
    presentGeneSequence = dict()
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
                print("\n".join([">" + k + "\n" + v + "\n" for k, v in presentGeneSequence.items()]), file=presentGeneName + ".multi.fasta")
                # We put the new informations in place
                presentGeneName = annotationGene
                presentgeneSens = annotationSens
                presentGeneSequence = dict()

        # If the annotationInterSca is after the scaffold of the vcf go to the next position of the vcf
        if (annotationIntervalSca > mutationScaffold):
            mutationPosition = increaseVCF(vcf_reader, annotationIntervalSca, annotationIntervalStart)
            mutationScaffold = mutationPosition.CHROM
            mutationPositionNumber = mutationPosition.POS
        # We have gone to the next interesting vcf position

        # If the bedInterSca is before the scaffold of the vcf go to the next position of the bed
        if (annotationIntervalSca < mutationScaffold):
            # ALL THE ANNOTATION HAVE NO CORRESPONDANT IN VCF, ENTIRE ANNOTATION WITH ONLY N FOR ALL BUT REFERENCE
            presentGeneSequence = addDic(presentGeneSequence, writeAbsentVCF(sampleNames, annotationIntervalStop - annotationIntervalStart))
            continue

        if (annotationIntervalStart < mutationPositionNumber):
            mutationPosition = increaseVCF(vcf_reader, annotationIntervalSca, annotationIntervalStart)
            mutationScaffold = mutationPosition.CHROM
            mutationPositionNumber = mutationPosition.POS
        # We have gone to the next interesting vcf position

        # The annotation interval is completely before the vcf position, we go to next group of position in annotation
        if (annotationIntervalStop < mutationPositionNumber):
            # ALL THE ANNOTATION HAVE NO CORRESPONDANT IN VCF, ENTIRE ANNOTATION WITH ONLY N FOR ALL BUT REFERENCE
            presentGeneSequence = addDic(presentGeneSequence, writeAbsentVCF(sampleNames, annotationIntervalStop - annotationIntervalStart))
            continue
        # The annotation interval is after the vcf position, we go to next position of vcf

        # The vcf position is in the BED interval
        if (mutationPositionNumber > annotationIntervalStart):
            presentGeneSequence = addDic(presentGeneSequence, writeAbsentVCF(sampleNames, mutationPositionNumber - annotationIntervalStart))
        while (mutationScaffold == annotationIntervalSca and mutationPositionNumber <= annotationIntervalStop):
            d = generatePosition(mutationPosition)
            addDic(presentGeneSequence, d)
            mutationPosition = increaseVCF(vcf_reader, annotationIntervalSca, annotationIntervalStart)
            mutationScaffold = mutationPosition.CHROM
            mutationPositionNumber = mutationPosition.POS


# BEGIN: OPENING
vcfHandle = open(args.input_vcf[0])
vcf_reader = vcf.Reader(vcfHandle)
annotationHandle = open(args.input_annotation3[0])

main()

# END: CLOSING
vcfHandle.close()
annotationHandle.close()


# for position in vcf_reader:
#    scaffold = position.CHROM
#    positionNumber = position.POS
