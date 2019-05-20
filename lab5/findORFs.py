#!/usr/bin/env python3
# Name: Eric Mockler (emockler)
# Group Members: None

########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'findORFs: A command-line tool to find open reading frames (ORFs) in DNA sequences', 
                                             epilog = 'Use multiple start or stop flags to append boundary codons.', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1 -option2 <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= range(99, 1001), default=99, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'], nargs='?', help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'], nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Execution
########################################################################
def main(inCL=None):
    '''
    Find some genes.  
    '''
    # load FastAreader
    import sequenceAnalysis
    from operator import itemgetter
    ORFreader = sequenceAnalysis.FastAreader()
    if inCL is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(inCL)
    # locate ORFs in sequences
    for head, seq in ORFreader.readFasta():
        foundORFs = sequenceAnalysis.ORFfinder(seq, myCommandLine.args.start, myCommandLine.args.stop)
        foundORFs.finder()
        foundORFs.ORFlist.sort(key = lambda x: (x[3], -x[1]), reverse = True)
        print(head)
        indexList = []
        for frame in foundORFs.ORFlist:
            if myCommandLine.args.longestGene:
                # find longest ORF in a presorted ORF list 
                if frame[0] < 0:
                    if frame[1] not in indexList and frame[3] >= myCommandLine.args.minGene:
                        indexList.append(frame[1])
                        print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame[0], frame[1], frame[2], frame[3]))
                else:
                    if frame[2] not in indexList and frame[3] >= myCommandLine.args.minGene:
                        indexList.append(frame[2])
                        print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame[0], frame[1], frame[2], frame[3]))                     
            else:
                if frame[3] >= myCommandLine.args.minGene:
                    print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame[0], frame[1], frame[2], frame[3]))
    
if __name__ == "__main__":
    main()  # delete the list when you want to run with STDIN


