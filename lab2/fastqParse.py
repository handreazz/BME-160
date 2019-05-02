#!/usr/bin/env python3 
# Name: Eric Mockler (emockler) 

'''
Parses FASTQ name information from a FASTQ-formatted line 

in: 
    @EAS139:136:FC706VJ:2:2104:15343:197393

out:
    Instrument = EAS139
    Run ID = 136
    Flow Cell ID = FC706VJ
    Flow Cell Lane = 2
    Tile Number = 2104
    X-coord = 15343
    Y-coord = 197393
'''

class FastqString(str):
    ''' Parses a FASTQ string.'''
    def length(self):
        ''' Returns length of input '''
        return len(self)
    def print(self):
        ''' Prints split FASTQ string by delimiting colons.'''   
        splitFASTQ = self.replace('@', '')
        splitFASTQ = splitFASTQ.replace(' ', '')
        splitFASTQ = splitFASTQ.split(':')
        print(
"""
Instrument = {0}
Run ID = {1}
Flow Cell ID = {2}
Flow Cell Lane = {3}
Tile Number  = {4}
X-coord = {5}
Y-coord = {6}
""".format(splitFASTQ[0], splitFASTQ[1], splitFASTQ[2], splitFASTQ[3],
            splitFASTQ[4], splitFASTQ[5], splitFASTQ[6])
        )

def main():
    ''' Print defined FASTQ name information.'''
    rawFASTQ = input("Input FASTQ line:")
    pureFASTQ = FastqString(rawFASTQ)
    pureFASTQ.print()
main()