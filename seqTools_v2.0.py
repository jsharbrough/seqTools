# -*- coding: utf-8 -*-
'''To run program in unix terminal, type the following:

	python seqTools_v2.0.py <command> <input>

in which command can be any of the following:

	buildCodonDict
			Produces aa sequences of nucleotide sequence input assuming the genetic code from user 
			input. Input should be in fasta format, output will aa sequences also in fasta format.
		Usage	
			python seqTools_v2.0.py buildCodonDict foo.fasta [options]
		Options
			-c	Dictates the genetic code to be used. [standard]
				vertebrateMt, echinodermMt, coelenterateMt, ascidianMt, standard, chlorophyceanMt, 
				yeastMt, bacterial, euplotidNuc, pterobranchiaMt, yeastNuc, flatwormMt, ciliateNuc, 
				invertebrateMt, trematodeMt

	buildSeqDict
			Primarily an internal tool that gets used to efficiently read-in fasta files. Will 
			return a key-value table with sequence tags as keys and sequences as values. 
		Usage
			python seqTools_v2.0.py buildSeqDict foo.fasta [options]
		Options

	reverseComplement
			Primarily an internal tool. If you would like to reverse translate all sequences in a 
			fasta file, use batchRevComp. reverseComplement can handle standard IUPAC ambiguity 
			codes (up to two alleles). Returns the reverse complement of the input sequence. 
		Usage
			python seqTools_v2.0.py reverseComplement <sequence> [options]
		Options

	cigarCalc
			A useful tool without much infrastructure. More development coming. This tool decodes 
			cigar strings produced by mapping programs and returns the length of the alignment. 
		Usage
			python seqTools_v2.0.py cigarCalc <cigar_string> [options]
		Options

	unwrap
			Tool to unwrap fasta file sequences. Takes a fasta file as input and removes 
			intra-sequence carriage returns and tabs.
		Usage
			python seqTools_v2.0.py unwrap foo.fasta [options]
		Options

	batchRevComp
			Tool that obtains the reverse complement sequence of all sequences in a file and 
			outputs them to a new file. Can handle all IUPAC ambiguity codes (up to two 
			nucleotides).
		Usage
			python seqTools_v2.0.py batchRevComp foo.fasta [options]
		Options

	batchTranslate
			Tool to translate many fasta files. Takes a file of file names as input and 
			sequenctially opens each file and translates it by calling the buildCodonDict 
			submodule.
		Usage
			python seqTools_v2.0.py batchTranslate foo.fofn [options]
		Options
	
	vcf2Fasta
			Tool that produces a fasta sequence from a Variant Call Formatted (VCF) file. Tool uses 
			standard IUPAC ambiguity codes (up to two alleles). Current filters include minimum 
			depth over the site and mapping quality assigned by the mapping program. 
		Usage
			python seqTools_v2.0.py vcf2Fasta foo.vcf [options]
		Options
			-rc	Reverse complement. Output will be the reverse complement of the mapped sequence.
				[False]
				True, False
			-d	Minimum depth requirement. Base calls with depth of coverage lower than -d will be 
				output as missing data (N). [10]
			-q	Minimum mapping quality threshold. Base calls with lower Phred-scaled mapping 
				quality than -q will be output as missing data (N). [30]

	help
			Writes this menu to the standard output.
			
Input should be in the form required by the submodule command. Note the program is case sensitive.

If you would like to know which command to use, type:

	python seqTools_v2.0.py help

Questions should be directed to Joel Sharbrough at jsharbro[at]gmail.com 

Copyright (c) 2017 Joel Sharbrough
'''
import sys
readMe = 'To run program in unix terminal, type the following:\n\n\tpython seqTools_v2.0.py <command> <input>\n\nin which command can be any of the following:\n\n\tbuildCodonDict\n\t\t\tProduces aa sequences of nucleotide sequence input assuming the genetic code from user \n\t\t\tinput. Input should be in fasta format, output will aa sequences also in fasta format.\n\t\tUsage\n\t\t\tpython seqTools_v2.0.py buildCodonDict foo.fasta [options]\n\t\tOptions\n\t\t\t-c\tDictates the genetic code to be used. [standard]\n\t\t\t\tvertebrateMt, echinodermMt, coelenterateMt, ascidianMt, standard, chlorophyceanMt, \n\t\t\t\tyeastMt, bacterial, euplotidNuc, pterobranchiaMt, yeastNuc, flatwormMt, ciliateNuc, \n\t\t\t\tinvertebrateMt, trematodeMt\n\n\tbuildSeqDict\n\t\t\tPrimarily an internal tool that gets used to efficiently read-in fasta files. Will \n\t\t\treturn a key-value table with sequence tags as keys and sequences as values. \n\t\tUsage\n\t\t\tpython seqTools_v2.0.py buildSeqDict foo.fasta [options]\n\t\tOptions\n\n\treverseComplement\n\t\t\tPrimarily an internal tool. If you would like to reverse translate all sequences in a \n\t\t\tfasta file, use batchRevComp. reverseComplement can handle standard IUPAC ambiguity \n\t\t\tcodes (up to two alleles). Returns the reverse complement of the input sequence. \n\t\tUsage\n\t\t\tpython seqTools_v2.0.py reverseComplement <sequence> [options]\n\t\tOptions\n\n\tcigarCalc\n\t\t\tA useful tool without much infrastructure. More development coming. This tool decodes \n\t\t\tcigar strings produced by mapping programs and returns the length of the alignment. \n\t\tUsage\n\t\t\tpython seqTools_v2.0.py cigarCalc <cigar_string> [options]\n\t\tOptions\n\n\tunwrap\n\t\t\tTool to unwrap fasta file sequences. Takes a fasta file as input and removes \n\t\t\tintra-sequence carriage returns and tabs.\n\t\tUsage\n\t\t\tpython seqTools_v2.0.py unwrap foo.fasta [options]\n\t\tOptions\n\n\tbatchRevComp\n\t\t\tTool that obtains the reverse complement sequence of all sequences in a file and \n\t\t\toutputs them to a new file. Can handle all IUPAC ambiguity codes (up to two \n\t\t\tnucleotides).\n\t\tUsage\n\t\t\tpython seqTools_v2.0.py batchRevComp foo.fasta [options]\n\t\tOptions\n\n\tbatchTranslate\n\t\t\tTool to translate many fasta files. Takes a file of file names as input and \n\t\t\tsequenctially opens each file and translates it by calling the buildCodonDict \n\t\t\tsubmodule.\n\t\tUsage\n\t\t\tpython seqTools_v2.0.py batchTranslate foo.fofn [options]\n\t\tOptions\n\n\tvcf2Fasta\n\t\t\tTool that produces a fasta sequence from a Variant Call Formatted (VCF) file. Tool uses \n\t\t\tstandard IUPAC ambiguity codes (up to two alleles). Current filters include minimum \n\t\t\tdepth over the site and mapping quality assigned by the mapping program. \n\t\tUsage\n\t\t\tpython seqTools_v2.0.py vcf2Fasta foo.vcf [options]\n\t\tOptions\n\t\t\t-rc\tReverse complement. Output will be the reverse complement of the mapped sequence.\n\t\t\t\t[False]\n\t\t\t\tTrue, False\n\t\t\t-d\tMinimum depth requirement. Base calls with depth of coverage lower than -d will be \n\t\t\t\toutput as missing data (N). [10]\n\t\t\t-q\tMinimum mapping quality threshold. Base calls with lower Phred-scaled mapping \n\t\t\t\tquality than -q will be output as missing data (N). [30]\n\n\thelp\n\t\t\tWrites this menu to the standard output.\n\t\t\t\nInput should be in the form required by the submodule command. Note the program is case sensitive.\n\nIf you would like to know which command to use, type:\n\n\tpython seqTools_v2.0.py help\n\nQuestions should be directed to Joel Sharbrough at jsharbro[at]gmail.com\n\nCopyright (c) 2017 Joel Sharbrough'


def buildCodonDict(fasta,code='standard'):
    #uses standard genetic code
    geneticCodes = {'standard':{"TTT":"F",	"TTC":"F",	"TTA":"L",	"TTG":"L",	"TCT":"S",	"TCC":"S",	"TCA":"S",	"TCG":"S",	"TAT":"Y",	"TAC":"Y",	"TAA":"*",	"TAG":"*",	"TGT":"C",	"TGC":"C",	"TGA":"*",	"TGG":"W",	"CTT":"L",	"CTC":"L",	"CTA":"L",	"CTG":"L",	"CCT":"P",	"CCC":"P",	"CCA":"P",	"CCG":"P",	"CAT":"H",	"CAC":"H",	"CAA":"Q",	"CAG":"Q",	"CGT":"R",	"CGC":"R",	"CGA":"R",	"CGG":"R",	"ATT":"I",	"ATC":"I",	"ATA":"I",	"ATG":"M",	"ACT":"T",	"ACC":"T",	"ACA":"T",	"ACG":"T",	"AAT":"N",	"AAC":"N",	"AAA":"K",	"AAG":"K",	"AGT":"S",	"AGC":"S",	"AGA":"R",	"AGG":"R",	"GTT":"V",	"GTC":"V",	"GTA":"V",	"GTG":"V",	"GCT":"A",	"GCC":"A",	"GCA":"A",	"GCG":"A",	"GAT":"D",	"GAC":"D",	"GAA":"E",	"GAG":"E",	"GGT":"G",	"GGC":"G",	"GGA":"G",	"GGG":"G"},'invertebrateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'vertebrateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': '*', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': '*', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'yeastMt':{'CTT': 'T', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'T', 'CTA': 'T', 'CTC': 'T', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'coelenterateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'ciliateNuc':{'CTT': 'L', 'TAG': 'Q', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': 'Q', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'echinodermMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'euplotidNuc':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'C', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'bacterial':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'yeastNuc':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'S', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'ascidianMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'G', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'G', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'flatwormMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': 'Y', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'chlorophyceanMt':{'CTT': 'L', 'TAG': 'L', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'trematodeMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'pterobranchiaMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'K', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}}
    geneticCode = geneticCodes[code]
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
    return scaffoldDict, scaffoldList
    
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
    seqDict,seqList = buildSeqDict(fasta)
    outfile = open(fasta[0:-6] + '_revComp.fasta','w')
    for seq in seqList:
        revSeq = reverseComplement(seqDict[seq])
        outfile.write(seq +'_complement\n' + revSeq + '\n')
    outfile.close()

def batchTranslate(fofn):
    infile = open(fofn,'r')
    for line in infile:
        fastaFile = line[0:-1]
        buildCodonDict(fastaFile)
    infile.close()

def vcf2Fasta(vcf,complement='False',minDepth='10',minQual='30'):
    if complement == 'True':
        complement = True
    else:
        complement = False
    minDepth = int(minDepth)
    minQual = int(minQual)
    infile = open(vcf,'r')
    outfile = open(vcf[0:-4] + '.fasta','w')
    outfile.write('>' + vcf[0:-4] + '\n')
    ambiguityCodes = {("A","C"):"M",("A","G"):"R",("A","T"):"W",("C","G"):"S",("C","T"):"Y",("G","T"):"K",("C","A"):"M",("G","A"):"R",("T","A"):"W",("G","C"):"S",("T","C"):"Y",("T","G"):"K"}
    seq = ''
    for line in infile:
        if line[0] != '#':
            lineSplit = line.split('\t')
            if lineSplit[5] != '.':
                qual = float(lineSplit[5])
            else:
                qual = 0
                print line
            chrom = lineSplit[0]
            ref = lineSplit[3]
            alt = lineSplit[4]
            info = lineSplit[7]
            gtInfo = lineSplit[9]
            infoSplit = info.split(';')
            gtInfoSplit = gtInfo.split(':')
            gt = gtInfoSplit[0]
            for item in infoSplit:
                if item[0:3] == 'DP=':
                    depth = int(item[3:])
            if gt == '0/0':
                nucleotide = ref
            elif gt == '1/1' and len(alt) == 1:
                nucleotide = alt
            elif gt == '0/1':    
                if len(ref) == len(alt):
                    nucleotide = ambiguityCodes[(ref,alt)]
                elif len(ref) >= len(alt):
                    nucleotide = ref
                else:
                    nucleotide = alt
            elif gt == '1/2':
                altSplit = alt.split(',')
                nucleotide = ambiguityCodes[(altSplit[0],altSplit[1])]
            else:
                nucleotide = 'N'
            if qual >= minQual and depth >= minDepth:
                seq += nucleotide
            else:
                seq += 'N'
    logfile = open('vcf2Fasta_logfile.txt','a')
    logfile.write(vcf[0:-4] + '\t' + str(len(seq)) + '\n')
    logfile.close()
    if complement == True:
        seq = reverseComplement(seq)
    outfile.write(seq + '\n')
    infile.close()
    outfile.close()

def help():
    sys.stdout.write(readMe)

def buildCommand():
    command = sys.argv[1]
    programList = ['buildCodonDict','buildSeqDict','reverseComplement','cigarCalc','unwrap','batchRevComp','batchTranslate','vcf2Fasta','help']
    programDict = {'batchRevComp': [], 'batchTranslate': [], 'buildCodonDict': ['-c'], 'buildSeqDict': [], 'cigarCalc': [], 'help': [], 'reverseComplement': [], 'unwrap': [], 'vcf2Fasta': ['-rc', '-d', '-q']}
    commandStatement = 'help()'
    if command == 'help' or command not in programDict:
        commandStatement = command + '()'
    else:
        Input = sys.argv[2]
        options = programDict[command]
        if len(options) == 0 or len(sys.argv) == 3:
            logfile = open('seqTools.log','a')
            logfile.write('Running ' + command + ' assuming default parameters\n')
            logfile.close()
            commandStatement = command + '("' + Input + '")'
        elif command == 'buildCodonDict':
            optionValues = ['standard']
            if sys.argv[-1] != optionValues[0]:
                optionValues[0] = sys.argv[-1]
            logfile = open('seqTools.log','a')
            logfile.write('Running ' + command + ' assuming the following parameters:\n\t')
            j = 0
            for item in options:
                logfile.write(item + ' ' + optionValues[j] + '\t')
                j += 1
            logfile.close() 
            commandStatement = command + '("' + Input + '","' + sys.argv[-1] + '")'
        elif command == 'vcf2Fasta':
            i = 0
            optionValues = [False,10,30]
            while i+4 < len(sys.argv):
                if sys.argv[i+3] == '-rc':
                    if sys.argv[i+4] == 'True':
                        optionValues[0] = True
                elif sys.argv[i+3] == '-d':
                    if sys.argv[i+4] != '10':
                        optionValues[1] = int(sys.argv[i+4])
                elif sys.argv[i+3] == '-q':
                    if sys.argv[i+4] != '30':
                        optionValues[2] = int(sys.argv[i+4])
                i += 1
            logfile = open('seqTools.log','a')
            logfile.write('Running ' + command + ' assuming the following parameters:\n\t')
            j = 0
            for item in options:
                logfile.write(item + ' ' + optionValues[j] + '\t')
                j += 1
            logfile.write('\n')
            logfile.close()
            commandStatement = command + '("' + Input + '","' + str(optionValues[0]) + '","' + str(optionValues[1]) + '","' + str(optionValues[2]) + '")'
    exec(commandStatement)
    
buildCommand()