#!/usr/bin/env python3
# Name: Eric Mockler (emockler)

import sequenceAnalysis
import math
import operator
def main ():
    # make sure to change this to use stdin
    myReader1 = sequenceAnalysis.FastAreader('testGenome.fa') 
    myNuc1 = sequenceAnalysis.NucParams()
    for head, seq in myReader1.readFasta():
        myNuc1.addSequence(seq)
    # make sure to change this to use stdin
    myReader2 = sequenceAnalysis.FastAreader('haloVolc1_1-genes.fa') 
    myNuc2 = sequenceAnalysis.NucParams()
    for head, seq in myReader2.readFasta():
        myNuc2.addSequence(seq)

    # find gc content
    nucComp1 = myNuc1.nucComposition()
    nucComp2 = myNuc2.nucComposition()
    nucCount1 = myNuc1.nucCount()
    nucCount2 = myNuc2.nucCount()
    GC1 = nucComp1.get('G') + nucComp2.get('C')
    GC2 = nucComp2.get('G') + nucComp2.get('C')
    log2 = math.log(((GC1/nucCount1) * 100) / \
                    ((GC2/nucCount2) * 100), 2)
    # print gc content
    print("GC content (log2(GC1/GC2)): {:.2f}\n".format(log2))
    # sort codons in alpha order, by Amino Acid
    sortedCodons = sorted(myNuc1.rnaCodonTable.items(), \
                    key = operator.itemgetter(1, 0))
    # calculate relative codon usage & frequency per codon and print
    print("codon : amino log2(usage1/usage2) (log2(count1/count2))")
    codonComp1 = myNuc1.codonComposition()
    codonComp2 = myNuc2.codonComposition()
    AAcomp1 = myNuc1.aaComposition()
    AAcomp2 = myNuc2.aaComposition()
    for codonAApair in sortedCodons:
        nuc = codonAApair[0]
        aa = codonAApair[1]
        relativeComp = math.log(
            (codonComp1.get(nuc)/ \
                AAcomp1.get(aa) / \
            (codonComp2.get(nuc)/ \
                AAcomp2.get(aa))), 2)
        count = math.log((codonComp1.get(nuc) / \
                codonComp2.get(nuc)), 2)
        print ('{:s}    :    {:s}     {:5.2f}          ({:6.2f})'.format(
            nuc, aa, relativeComp*100, count))
    print('\r')
if __name__ == "__main__":
    main()
    