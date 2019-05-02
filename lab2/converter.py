#!/usr/bin/env python3 
# Name: Eric Mockler (emockler)

'''
Converts amino acid representations between the 3-letter DNA/RNA codon code, 
3-letter amino acid code, and 1-letter amino acid code.

in:             in:             in:
    atg             e               asp
out:            out:            out:
    ATG = MET       E = GLU         ASP = D

'''
short_AA = {
            'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
            }

long_AA = {value:key for key,value in short_AA.items()}

RNA_codon_table = {
# Second Base
# U             C             A             G
#U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
#C 
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
#A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
#G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
}

dnaCodonTable = {key.replace('U','T'):value for key, value in RNA_codon_table.items()}

nucleotides = ('G', 'T', 'A', 'C', 'U')

class AminoFormatConverter(str):
    """ Parses user-inputted amino code"""
    def length(self):
        """Returns length of input code"""
        return len(self)
    def parseInput(self):
        """ Determines type of user input and prints its complement output """
        aminoRaw = self.upper()
        if self.length() == 3:
            # if input is 3-letter amino
            if (aminoRaw[0] not in nucleotides) and (aminoRaw[1] not in nucleotides) and (aminoRaw[2] not in nucleotides):
                self.shortAAGet(aminoRaw)
            # if input is 3-letter amino that contains a conflicting nucleotide abbreviation in 2nd char position
            elif (aminoRaw[0] not in nucleotides) or (aminoRaw[2] not in nucleotides) and (aminoRaw[1] in nucleotides):
                self.shortAAGet(aminoRaw)
            # if input is 3-nucleotide codon code    
            elif (aminoRaw[0] in nucleotides) and (aminoRaw[1] in nucleotides) and (aminoRaw[2] in nucleotides):
                if 'U' in aminoRaw:
                    self.longRNAGet(aminoRaw)
                elif 'T' in aminoRaw:
                    self.longDNAGet(aminoRaw)
        elif self.length() == 1:
            self.longAAGet(aminoRaw)
        else:
            print("{0} = unknown".format(aminoRaw))

    def shortAAGet(self, amino):
        """ Prints short AA code from long code"""
        shortCode = short_AA.get(amino, "unknown")
        print("{0} = {1}".format(amino, shortCode))
    def longAAGet(self, amino):
        """ Prints long AA code from short AA code """
        longCode = long_AA.get(amino, "unknown")
        print("{0} = {1}".format(amino, longCode))
    def longRNAGet(self, amino):
        """ Prints long AA code from RNA codon code """
        shortCode = RNA_codon_table.get(amino, "unknown")
        print("{0} = {1}".format(amino, shortCode))
    def longDNAGet(self, amino):
        """ Prints long AA code from DNA codon code """
        shortCode = dnaCodonTable.get(amino, "unknown")
        print("{0} = {1}".format(amino, shortCode))
    
def main():
    ''' Instructs user input of amino acid code'''
    rawAmino = input("Enter a DNA/RNA codon code, 3-letter amino acid, or 1-letter amino acid: ")
    rawAmino = AminoFormatConverter(rawAmino)
    rawAmino.parseInput()
main()