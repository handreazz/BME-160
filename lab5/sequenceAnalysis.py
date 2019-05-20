#!/usr/bin/env python3
# Name: Eric Mockler (emockler)

class ProteinParam :
    """ 
    Return parameters of an input protein sequence. 

    Instanstiation: 
        proteinParams = ProteinParam(protein_sequence)
    Usage: 
        proteinPI = proteinParams.pI([floating_point_precision])
    """

# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        """ Initialize AA composition dict of input string """
        self.rawSeq = protein.upper()
        self.aa2mwSeq = {key: 0 for key in ProteinParam.aa2mw}
        self.aa2repSeq = {key: 0 for key in ProteinParam.aa2mw}
        self.cleanSeq = ''
        # update aa2mwSeq with molar weight if amino acid exists in protein
        for c in self.rawSeq:
            if c in ProteinParam.aa2mw:
                self.aa2mwSeq.update({c: ProteinParam.aa2mw.get(c)})
                self.cleanSeq += c
            # count repeated AA
            if self.aa2mwSeq.get(c):
                self.aa2repSeq.update({c: (self.aa2repSeq.get(c) + 1)})
     
    def _charge_ (self, pH):
        """ Return the net charge of the input protein sequence with respect to pH """
        sumChargePos = sumChargeNeg = 0
        for c in ProteinParam.aa2chargePos:
            if (self.aa2repSeq.get(c)):
                pKaPos = ProteinParam.aa2chargePos.get(c)
                sumChargePos += self.aa2repSeq.get(c) * ((10 ** pKaPos) / ((10 ** pKaPos) + (10 ** pH))) 
        # Add N-terminus to positive charge
        sumChargePos += (10 ** ProteinParam.aaNterm) / ((10 ** ProteinParam.aaNterm) + (10 ** pH))
        for c in ProteinParam.aa2chargeNeg:
            if (self.aa2repSeq.get(c)):
                pKaNeg = ProteinParam.aa2chargeNeg.get(c)
                sumChargeNeg += self.aa2repSeq.get(c) * ((10 ** pH) / ((10 ** pKaNeg) + (10 ** pH)))
        # Add C-terminus to positive charge
        sumChargeNeg += (10 ** pH) / ((10 ** ProteinParam.aaCterm) + (10 ** pH))
        return sumChargePos - sumChargeNeg
    
    def _setPrecision_ (self, precision):
        """ Return an iterator reflecting the desired floating point precision. """
        decimalPlace = '0.'
        for i in range(precision):  
            if i == precision - 1: decimalPlace += '1'
            else: decimalPlace += '0'
            i += 1        
        return float(decimalPlace)

    def _binarySearch_ (self, pHRange, p, r, target):
        """ Binary search for a pH keyed to a target charge of ~0 """
        q = (p+r) // 2
        pHChargePair = pHRange[q]
        charge = pHChargePair[1]
        pH = pHChargePair[0]
        # recursive base case
        if p > r:
            return pH
        # narrow search space to left or right
        else:
            if target > charge >= 0:
                return pH
            elif target < charge:
                return self._binarySearch_(pHRange, p, q-1, target)
            elif target > charge:
                return self._binarySearch_(pHRange, q+1, r, target)
    # end of private methods

    def aaCount (self):
        """ Sum the number of AA in input protein """
        return len(self.cleanSeq)

    def pI (self, precision = 2):
        """ Returns the theoretical isolelectric point of the input protein sequence """
        pH = 0.0
        pHRange = []
        preciseIterator = self._setPrecision_(precision)
        while pH <= 14.0:
            pHRange.append([pH, 0]) 
            pH += preciseIterator
        for i, pHChargePair in enumerate(pHRange):
            charge = self._charge_(pHChargePair[0])
            pHRange[i] = [pHChargePair[0], charge]
        # Binary search for pH with charge between preciseIterator and 0
        pHRange = sorted(pHRange, key = lambda charge:charge[1])
        pHRangeLen = len(pHRange)
        pI = self._binarySearch_(pHRange, 0, pHRangeLen, preciseIterator)  
        return pI   

    def aaComposition (self) :
        """ Return dictionary containing composition of AA found in input protein """
        return self.aa2repSeq

    def molarExtinction (self, cystine = True):
        """ Return molar extinction coefficient of input protein """
        if cystine:
            return (self.aa2repSeq.get('Y') * ProteinParam.aa2abs280.get('Y')) + (self.aa2repSeq.get('W') 
                    * ProteinParam.aa2abs280.get('W')) + (self.aa2repSeq.get('C') * ProteinParam.aa2abs280.get('C'))
        else:
            return (self.aa2repSeq.get('Y') * ProteinParam.aa2abs280.get('Y')) + (self.aa2repSeq.get('W') 
                    * ProteinParam.aa2abs280.get('W'))

    def massExtinction (self, cystine = True):
        """ Calculate the mass extinction coefficient of the input AA sequence """
        myMW =  self.molecularWeight()
        return self.molarExtinction(cystine) / myMW if myMW else 0.0

    def molecularWeight (self):
        """ Calculate the molecular weight of the input protein """
        sumMass = 0
        for c in self.aa2mwSeq:
            if self.aa2mwSeq.get(c):
                sumMass += self.aa2repSeq.get(c) * (self.aa2mwSeq.get(c) - ProteinParam.mwH2O)
        return ProteinParam.mwH2O + sumMass

class NucParams :
    """
    Track codon counts in a genetic sequence
        initialization:
            protein = NucParams(DNA_or_RNA_sequence)
        usage:
            proteinCodonComp = protein.codonComposition() 
    """
            
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        """ 
        constructor: store set of nucleotides for comparison 
        & initialize codon dictionary from input sequence 
        """
        self.nucSet = {'A', 'C', 'T', 'G', 'U', 'N'}
        self.aaComp = {amino: 0 for amino in set(NucParams.rnaCodonTable.values())}
        self.nucComp = {nuc: 0 for nuc in self.nucSet}
        self.codonComp = {codon: 0 for codon in NucParams.rnaCodonTable.keys()}
        self.aaList = []
        self.codonList = []
        self.addSequence(inString)

    def addSequence (self, inSeq):
        """ append codon sequence to codon dictionary """
        # if inSeq is not empty
        if inSeq:
            for nuc in inSeq:
                if nuc not in self.nucSet:
                    inSeq.replace(nuc, '')
            inSeqCodonList = [inSeq[j:j+3] for j in range(0, len(inSeq), 3)]
            self.codonList.extend(inSeqCodonList)
            # translate AA sequence from codons
            for codon in inSeqCodonList:
                if 'U' in codon:
                    self.aaList.append(self.rnaCodonTable.get(codon))
                else:
                    self.aaList.append(self.dnaCodonTable.get(codon))
   
    def aaComposition(self):
        """ return AA composition of codon dictionary """
        for amino in self.aaList:
            self.aaComp.update({amino: (self.aaComp.get(amino) + 1)})
        return self.aaComp
   
    def nucComposition(self):
        """ return composition of valid nucleotides in dictionary """
        for codon in self.codonList:
            for nuc in codon:
                self.nucComp.update({nuc: (self.nucComp.get(nuc) + 1)})
        return self.nucComp
   
    def codonComposition(self): 
        """ return composition of found codons in dictionary """
        for codon in self.codonList:    
           rnaCodon = codon.replace('T', 'U')
           self.codonComp.update({rnaCodon: (self.codonComp.get(rnaCodon) + 1)})
        return self.codonComp
   
    def nucCount(self):
        """ return count of valid nucleotides """
        return len(self.aaList)*3

class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
    print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            import sys
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()

            while not line.startswith('>') :
                if not line: # we are at EOF
                    return header, sequence
                line = fileH.readline()

            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

class ORFfinder:
    """
    Find open reading frames (ORFs) in a double-stranded DNA sequence, in all 
    6 reading frames on the sense and antisense strands. Start and/or 
    stop codons denote an ORF, and are defined by the user.

    initialization:
    ORFs = ORFfinder(DNA sequence[, start codon list[, stop codon list]])
    usage:
    print(ORFs.ORFlist[i]) -> (frame, start, stop, length)
    """
    def __init__(self, nucSequence,
        start = ['ATG'], stop = ['TAG','TGA','TAA']):
        """ 
        constructor: saves DNA sequence 
        & reference start/stop codons
        """
        self.nucSequence = nucSequence
        self.complementSeq = self._complement_(self.nucSequence)
        self.startCodons = start
        self.stopCodons = stop
        self.ORFlist = []
        # lookup table for preexisting starts or stops 
        self.foundStartorStop = []
    
    def _complement_(self, seq):
        """ return the complement of a DNA sequence (without reversing stream) """
        complementSeq = ''
        for c in seq:
            if c is 'A': complementSeq += 'T'
            elif c is 'T': complementSeq += 'A'
            elif c is 'C': complementSeq += 'G'
            elif c is 'G': complementSeq += 'C'
        return complementSeq

    def _reverseIndex_(self, index, indexParity):
        """ return index that starts iteration at the end of the complement
        sequence, reversing stream to enable iteration over antisense strand """
        if indexParity is '-':
            return index % len(self.nucSequence)
        else: 
            return -index - 1
    # end of private methods  

    def finder(self):
        """ :
        find ORFs by iterating through each potential
        reading frame
        """
        # set found start and stop codon flags for sense & antisense
        self.foundStartorStop = [[frame, [False, False], [False, False]] for frame in range(1, 4)]
        # find ORFs in frames 1, 2, 3 on sense & antisense strands
        for i in range(0, len(self.nucSequence), 3):
            self.findInFrame(i, 1)
            self.findInFrame(i, 2)
            self.findInFrame(i, 3)
        # find frames with no start/stop codons to reflect in ORFList
        for foundFrame in self.foundStartorStop:
            if True not in foundFrame[1]:
                self.ORFlist.append((foundFrame[0], 1, len(self.nucSequence), len(self.nucSequence)))
            if True not in foundFrame[2]:
                 self.ORFlist.append((foundFrame[0] * -1, 1, len(self.nucSequence), len(self.nucSequence)))

    def findInFrame(self, index, frame):
        """
        check if a codon at a specified index & frame match start or stop
        codons, denoting an ORF to be reported
        """
        framedIndex = index + frame - 1 
        currentCodon = self.nucSequence[framedIndex:framedIndex+3]
        if len(currentCodon) < 3: return
        if currentCodon in self.startCodons:
            self.foundStartorStop[frame-1][1][0] = True
            self.geneReporter(framedIndex, frame, True)
        elif currentCodon in self.stopCodons and framedIndex + 3 < len(self.nucSequence):
            self.foundStartorStop[frame-1][1][1] = True
            self.geneReporter(framedIndex, frame, False)
        # check reverse complement
        frame = -frame
        reverseIndex = self._reverseIndex_(framedIndex, '+')
        if reverseIndex is -1:
            currentCodon = self.complementSeq[reverseIndex-2:] 
        else:
            currentCodon = self.complementSeq[reverseIndex-2:reverseIndex+1]
        currentCodon = currentCodon[::-1]
        if currentCodon in self.startCodons:
            foundORFreverse = True
            self.foundStartorStop[-frame-1][2][0] = True  
            self.geneReporter(reverseIndex, frame, True, foundORFreverse)
        elif currentCodon in self.stopCodons and reverseIndex - 3 > -len(self.nucSequence):
            foundORFreverse = True           
            self.foundStartorStop[-frame-1][2][1] = True 
            self.geneReporter(reverseIndex, frame, False, foundORFreverse)
    
    def geneReporter(self, framedIndex, frame,
    foundStartCodon, foundORFreverse = False):
        """ report all ORFs in the stream starting at a start or stop codon """
        seq = self.complementSeq if foundORFreverse else self.nucSequence
        boundaryCodons = self.stopCodons if foundStartCodon else self.startCodons
        previousStartCodon = self.foundStartorStop[-frame-1][2][0] if foundORFreverse else self.foundStartorStop[frame-1][1][0]
        if foundORFreverse: 
            frame = -frame
            framedIndex = self._reverseIndex_(framedIndex, '-')
            boundary = 0
            iterator = -3
        else:
            boundary = len(self.nucSequence) - 1
            iterator = 3
        start = stop = framedIndex        
        # if start found, move framedIndex downstream until a stop codon is found
        if foundStartCodon is True:
            for i in range(framedIndex, boundary, iterator):
                if foundORFreverse: # if complement strand, reverse antisense codon
                    currentCodon = seq[i-2:] if i is -1 else seq[i-2:i+1]
                    currentCodon = currentCodon[::-1]
                else: 
                    currentCodon = seq[i:i+3]
                stop = i
                if currentCodon in boundaryCodons:
                    start = start + 1
                    if foundORFreverse:
                        self.foundStartorStop[frame-1][2][1] = True
                        frame = -frame
                        stop = stop - 1
                        self.ORFlist.append((frame, stop, start, (start - stop) + 1))
                        break
                    else:
                        stop = stop + 3
                        self.foundStartorStop[frame-1][1][1] = True
                        self.ORFlist.append((frame, start, stop, (stop - start) + 1))
                        break
            # if no stop codon found, report boundary case
            if foundORFreverse: boundaryExists = self.foundStartorStop[-frame-1][2][1]  
            else: boundaryExists = self.foundStartorStop[frame-1][1][1]
            if not boundaryExists:
                start = start + 1
                stop = 1 if foundORFreverse else len(self.nucSequence)
                if foundORFreverse: 
                    frame = -frame
                    self.ORFlist.append((frame, stop, start, start - stop))
                else:
                    self.ORFlist.append((frame, start, stop, stop - start))
        # if stop found and start is left hanging, report boundary case  
        elif not previousStartCodon:
            if foundORFreverse: 
                frame = -frame
                start = start - 1
                stop = len(self.nucSequence)
            else:
                stop = stop + 3
                start = 1
            self.ORFlist.append((frame, start, stop, (stop - start) + 1))
   