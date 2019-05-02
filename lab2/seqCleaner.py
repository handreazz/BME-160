#!/usr/bin/env python3 
# Name: Eric Mockler (emockler) 

'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example: 
 input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''

class DNAstring (str):
    def length (self):
        return (len(self))
    
    def purify(self):
        ''' Return an upcased version of the string, collapsing a single run of Ns.'''
        upperDNA = self
        upperDNA = upperDNA.upper()
        
        """ For each character, determine if N sequence exists,
            then slice out and replace with N-count """
        firstNFound = False
        stopIndex = 0
        startIndex = 0
        for i, c in enumerate(upperDNA):
            if not c == 'N':
                firstNFound = False
            if not firstNFound and c == 'N':
                firstNFound = True
                startIndex = i
                if stopIndex == 0:
                    stopIndex = startIndex
            if c == 'N':
                stopIndex += 1
        countN = (stopIndex - startIndex)
        """ Assemble purified string """
        purifiedDNA = upperDNA[:startIndex] + "{" + str(countN) + "}" + upperDNA[stopIndex:]            
        return purifiedDNA
            
def main():
     
    ''' Get user DNA data and clean it up.'''
    data = input('DNA data?')
    thisDNA = DNAstring (data)
    pureData = thisDNA.purify()
    print (pureData)
    
main()