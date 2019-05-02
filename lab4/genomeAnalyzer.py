#!/usr/bin/env python3
# Name: Eric Mockler (emockler)

import sequenceAnalysis
import operator
def main ():
    # make sure to change this to use stdin
    myReader = sequenceAnalysis.FastAreader()
    myNuc = sequenceAnalysis.NucParams()
    for head, seq in myReader.readFasta():
        myNuc.addSequence(seq)
  
    # print sequence length
    print("\nsequence length: {:.2f} Mb\n".format(
        myNuc.nucCount() / (10**6)))
    # find gc content
    nucComp = myNuc.nucComposition()
    nucCount = myNuc.nucCount()
    GC = nucComp.get('G') + nucComp.get('C')
    # print gc content
    print("GC content: {:.1f}%\n".format((GC/nucCount) * 100))
    # sort codons in alpha order, by Amino Acid
    sortedCodons = sorted(myNuc.rnaCodonTable.items(), \
                    key = operator.itemgetter(1, 0))
    # calculate relative codon usage for each codon and print
    codonComp = myNuc.codonComposition()
    AAcomp = myNuc.aaComposition()
    for codonAApair in sortedCodons:
        nuc = codonAApair[0]
        aa = codonAApair[1]
        val = codonComp.get(nuc)/AAcomp.get(aa)
        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(
             nuc, aa, val*100, codonComp.get(nuc)))
    print('\r')
if __name__ == "__main__":
    main()
    