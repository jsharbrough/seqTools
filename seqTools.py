# -*- coding: utf-8 -*-
import sys
readMe = 'To run program in unix terminal, type the following:\n\n\tpython seqTools.py <command> <input>\n\nin which command can be any of the following:\n\n\tbuildCodonDict --Produces aa sequences of nucleotide sequence input using the standard genetic code. Input should be in fasta format, output will aa sequences also in fasta format. Execute by typing: python buildCodonDict foo.fasta\n\n\tbuildSeqDict -- Primarily an internal tool that gets used to efficiently read-in fasta files. Will return a key-value table with sequence tags as keys and sequences as values. Execute by typing: python buildSeqDict foo.fasta\n\n\treverseComplement --Primarily an internal tool. If you would like to reverse translate all sequences in a fasta file, use batchRevComp. reverseComplement can handle standard IUPAC ambiguity codes (up to two alleles). Returns the reverse complement of the input sequence. Execute by typing python reverseComplement <sequence>\n\n\tcigarCalc --A useful tool without much infrastructure. More development coming. This tool decodes cigar strings produced by mapping programs and returns the length of the alignment. Execute by typing python cigarCalc <cigar> where <cigar> is the cigar string of interest.\n\n\tunwrap --Tool to unwrap fasta file sequences. Takes a fasta file as input and removes intra-sequence carriage returns and tabs. Execute by typing python unwrap foo.fasta\n\n\tbatchRevComp --Tool that obtains the reverse complement sequence of all sequences in a file and outputs them to a new file. Can handle all IUPAC ambiguity codes (up to two nucleotides). Execute by typing: python batchRevComp foo.fasta\n\n\tbatchTranslate --Tool to translate many fasta files. Takes a file of file names as input and sequenctially opens each file and translates it by calling the buildCodonDict submodule. Execute by typing: python batchTranslate foo.fofn\n\n\thelp --Writes this menu to the standard output.\nand input should be in the form required by the submodule command.\n\nIf you would like to know which command to use, type:\n\n\tpython seqTools.py help'

programList = ['buildCodonDict','buildSeqDict','reverseComplement','cigarCalc','unwrap','batchRevComp','batchTranslate','help']


def buildCodonDict(fasta):
    #uses standard genetic code
    geneticCode = {'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': '*', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': '*', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}
    startCodons = ['ATT','ATC','ATA','ATG','GTG']
    seqDict = buildSeqDict(fasta)
    codonDict = {}
    AADict = {}
    for seq in seqDict:
        nucleotideSeq = seqDict[seq]
        codonList = []
        i = 2
        while i < len(nucleotideSeq):
            currCodon = nucleotideSeq[i-2] + nucleotideSeq[i-1] + nucleotideSeq[i]
            codonList.append(currCodon)
            i += 3
        codonDict[seq] = codonList
        AAseq = ''
        codonNum = 1
        for codon in codonList:
            if codonNum == 1:
                if codon in startCodons:
                    aa = 'M'
                else:
                    aa = geneticCode[codon]
            elif codon in geneticCode:
                aa = geneticCode[codon]
            else:
                aa = 'X'
            AAseq += aa
            codonNum += 1
        if AAseq[-1] == '*':
            AAseq = AAseq[0:-1]
        AADict[seq] = AAseq
    outfile = open(fasta[0:-6] + '_aaSeqs.fasta','w')
    for seq in seqDict:
        outfile.write(seq + '\n' + AADict[seq] + '\n')
    outfile.close()
    
def buildSeqDict(fasta):
    infile = open(fasta,'r')
    scaffoldDict = {}
    scaffoldList = []
    seqName = ''
    currSeq = ''
    for line in infile:
        if line[0] == '>':
            if seqName != '':
                scaffoldDict[seqName] = currSeq
            seqName = line
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            scaffoldList.append(seqName)
            currSeq = ''
        else:
            currSeq += line
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    scaffoldDict[seqName] = currSeq 
    return scaffoldDict
    
def reverseComplement(seq):
    #seqDict = buildSeqDict(fasta)
    #for sequence in seqDict:
        #seq = seqDict[sequence]
    seq_revc = ''
    for nuc in seq:
        if nuc == 'A':
            seq_revc = 'T' + seq_revc
        elif nuc == 'T':
            seq_revc = 'A' + seq_revc
        elif nuc == 'C':
            seq_revc = 'G' + seq_revc
        elif nuc == 'G':
            seq_revc = 'C' + seq_revc
        elif nuc == 'M':
            seq_revc = 'K' + seq_revc
        elif nuc == 'R':
            seq_revc = 'Y' + seq_revc
        elif nuc == 'S':
            seq_revc = 'S' + seq_revc
        elif nuc == 'W':
            seq_revc = 'W' + seq_revc
        elif nuc == 'K':
            seq_revc = 'M' + seq_revc
        elif nuc == 'Y':
            seq_revc = 'R' + seq_revc
        else:
            seq_revc = nuc + seq_revc
    return seq_revc

def cigarCalc(cigar):
    totalLength = 0
    currString = ''
    for char in cigar:
        if char == 'M':
            totalLength += int(currString)
            currString = ''
        elif char == 'N':
            totalLength += int(currString)
            currString = ''
        elif char == 'I':
            totalLength += int(currString)
            currString = ''
        elif char == 'D':
            currString = ''
        elif char == 'S':
            totalLength += int(currString)
            currString = ''
        elif char == 'H':
            currString = ''
        elif char == 'P':
            currString = ''
        elif char == 'X':
            totalLength += int(currString)
            currString = ''
        elif char == '=':
            totalLength += int(currString)
            currString = ''
        else:
            currString += char
    return totalLength

def unwrap(fasta):
    infile = open(fasta,'r')
    outfile = open(fasta[0:-6]+'_unwrapped.fasta','w')
    lineNum = 0
    for line in infile:
        if line[0] == '>':
            if lineNum == 0:
                outfile.write(line)
                lineNum += 1
            else:
                outfile.write('\n' + line)
        else:
            outfile.write(line[0:-1])
    infile.close()
    outfile.close()
    
def batchRevComp(fasta):
    fastaFile = open(fasta,'r')
    revSeqName, revSeq = reverseComplement(fastaFile)
    outfile = open(fastaFile[0:-6] + '_revComp.fasta','w')
    outfile.write(revSeqName +'\n' + revSeq + '\n')
    outfile.close()

def batchTranslate(fofn):
    infile = open(fofn,'r')
    for line in infile:
        fastaFile = line[0:-1]
        buildCodonDict(fastaFile)
    infile.close()

def help():
    sys.stdout.write(readMe)


if sys.argv[1] == 'help':
    help()
elif sys.argv[1] in programList:
    sys.argv[1](sys.argv[2])
else:
    help()