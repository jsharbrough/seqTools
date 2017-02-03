# -*- coding: utf-8 -*-
'''To run program in unix terminal, type the following:

	python seqTools_v2.3.py <command> <input>

in which command can be any of the following:

	buildCodonDict
			Produces aa sequences of nucleotide sequence input assuming the genetic code from user 
			input. Input should be in fasta format, output will aa sequences also in fasta format.
		Usage	
			python seqTools_v2.3.py buildCodonDict foo.fasta [options]
		Options
			-c	Dictates the genetic code to be used. [standard]
				vertebrateMt, echinodermMt, coelenterateMt, ascidianMt, standard, chlorophyceanMt, 
				yeastMt, bacterial, euplotidNuc, pterobranchiaMt, yeastNuc, flatwormMt, ciliateNuc, 
				invertebrateMt, trematodeMt

	buildSeqDict
			Primarily an internal tool that gets used to efficiently read-in fasta files. Will 
			return a key-value table with sequence tags as keys and sequences as values. 
		Usage
			python seqTools_v2.3.py buildSeqDict foo.fasta [options]
		Options

	reverseComplement
			Primarily an internal tool. If you would like to reverse translate all sequences in a 
			fasta file, use batchRevComp. reverseComplement can handle standard IUPAC ambiguity 
			codes (up to two alleles). Returns the reverse complement of the input sequence. 
		Usage
			python seqTools_v2.3.py reverseComplement <sequence> [options]
		Options

	cigarCalc
			A useful tool without much infrastructure. More development coming. This tool decodes 
			cigar strings produced by mapping programs and returns the length of the alignment. 
		Usage
			python seqTools_v2.3.py cigarCalc <cigar_string> [options]
		Options

	unwrap
			Tool to unwrap fasta file sequences. Takes a fasta file as input and removes 
			intra-sequence carriage returns and tabs.
		Usage
			python seqTools_v2.3.py unwrap foo.fasta [options]
		Options

	batchRevComp
			Tool that obtains the reverse complement sequence of all sequences in a file and 
			outputs them to a new file. Can handle all IUPAC ambiguity codes (up to two 
			nucleotides).
		Usage
			python seqTools_v2.3.py batchRevComp foo.fasta [options]
		Options

	batchTranslate
			Tool to translate many fasta files. Takes a file of file names as input and 
			sequenctially opens each file and translates it by calling the buildCodonDict 
			submodule.
		Usage
			python seqTools_v2.3.py batchTranslate foo.fofn [options]
		Options
	
	vcf2Fasta
			Tool that produces a fasta sequence from a Variant Call Formatted (VCF) file. Tool uses 
			standard IUPAC ambiguity codes (up to two alleles). Current filters include minimum 
			depth over the site and mapping quality assigned by the mapping program. 
		Usage
			python seqTools_v2.3.py vcf2Fasta foo.vcf [options]
		Options
			-rc	Reverse complement. Output will be the reverse complement of the mapped sequence.
				[False]
				True, False
			-d	Minimum depth requirement. Base calls with depth of coverage lower than -d will be 
				output as missing data (N). [10]
			-q	Minimum mapping quality threshold. Base calls with lower Phred-scaled mapping 
				quality than -q will be output as missing data (N). [30]
	
        seqRenamer
                        Tool to rename sequences within a fasta file. Input must be fasta. Output will be another 
                        fasta file with "renamed.fasta" as its suffix, and the original file name as its prefix.
                        Renaming scheme can be numeric (default) or a fixed name. Numeric renaming will add "_n"
                        to the end of the sequence tag, where n = the order of the sequence in the input fasta. 
                        Fixed renaming will add a user-specified string to the beginning of each sequence tag 
                        like follows: "><STRING>_sequenceTag".
                Usage
                        python seqTools_v2.3.py seqRenamer foo.fasta [options]
                Options
                        -n      Renaming method. If the default numeric renaming system is not used, the value
                                following the option flag will be appended as a prefix to all sequences in the 
                                fasta file. [numeric]
        
        fasta2fastq
                        Makes a fastq file with fake quality scores (ASCII_BASE 33) of -q from a fasta input. 
                        Note: input file handle MUST BE ".fasta"
                Usage
                        python seqTools_v2.3.py fastq2Fastq foo.fasta [options]
                Options
                        -q      Phred quality score for fake quality scores. Must be an integer. [40]
        
        nexus2Fasta
                        Tool for converting nexus files into fasta files. Input should be a nexus file. Output 
                        will be a fasta file.
                Usage
                        python seqTools_v2.3.py nexus2Fasta foo.nexus [options]
                Options
        
        fasta2Nexus
                        Tool for converting fasta files to nexus format. Input should be a fasta file, output 
                        will be a nexus file. Nexus file will be compatible with Paup or MrBayes, depending on 
                        option selected. Additional compatibility can be added by contacting the developer.
                Usage
                        python seqTools_v2.3.py fasta2Nexus foo.fasta [options]
                Options
                        -p      Identifies the compatibility mode desired. If you'd like additional compatibility 
                                added, please email developer. If paup compatibility is specified, a paup block at 
                                the bottom of the nexus file will be written, which instructs paup to construct a 
                                Neighbor Joinging Tree and a ML pairwise distance matrix of the alignment. [none]
                                none,mrbayes, paup
                                
        help
			Writes this menu to the standard output.
		Usage
		        python seqTools_v2.3.py help
			
Input should be in the form required by the submodule command. Note the program is case sensitive.

If you would like to know which command to use, type:

	python seqTools_v2.3.py help

Questions should be directed to Joel Sharbrough at jsharbro[at]gmail.com 

Copyright (c) 2017 Joel Sharbrough
'''
import sys
readMe = 'To run program in unix terminal, type the following:\n\n\tpython seqTools_v2.3.py <command> <input>\n\nin which command can be any of the following:\n\n\tbuildCodonDict\n\t\t\tProduces aa sequences of nucleotide sequence input assuming the genetic code from user \n\t\t\tinput. Input should be in fasta format, output will aa sequences also in fasta format.\n\t\tUsage\n\t\t\tpython seqTools_v2.3.py buildCodonDict foo.fasta [options]\n\t\tOptions\n\t\t\t-c\tDictates the genetic code to be used. [standard]\n\t\t\t\tvertebrateMt, echinodermMt, coelenterateMt, ascidianMt, standard, chlorophyceanMt, \n\t\t\t\tyeastMt, bacterial, euplotidNuc, pterobranchiaMt, yeastNuc, flatwormMt, ciliateNuc, \n\t\t\t\tinvertebrateMt, trematodeMt\n\n\tbuildSeqDict\n\t\t\tPrimarily an internal tool that gets used to efficiently read-in fasta files. Will \n\t\t\treturn a key-value table with sequence tags as keys and sequences as values. \n\t\tUsage\n\t\t\tpython seqTools_v2.3.py buildSeqDict foo.fasta [options]\n\t\tOptions\n\n\treverseComplement\n\t\t\tPrimarily an internal tool. If you would like to reverse translate all sequences in a \n\t\t\tfasta file, use batchRevComp. reverseComplement can handle standard IUPAC ambiguity \n\t\t\tcodes. Returns the reverse complement of the input sequence. \n\t\tUsage\n\t\t\tpython seqTools_v2.3.py reverseComplement <sequence> [options]\n\t\tOptions\n\n\tcigarCalc\n\t\t\tA useful tool without much infrastructure. More development coming. This tool decodes \n\t\t\tcigar strings produced by mapping programs and returns the length of the alignment. \n\t\tUsage\n\t\t\tpython seqTools_v2.3.py cigarCalc <cigar_string> [options]\n\t\tOptions\n\n\tunwrap\n\t\t\tTool to unwrap fasta file sequences. Takes a fasta file as input and removes \n\t\t\tintra-sequence carriage returns and tabs.\n\t\tUsage\n\t\t\tpython seqTools_v2.3.py unwrap foo.fasta [options]\n\t\tOptions\n\n\tbatchRevComp\n\t\t\tTool that obtains the reverse complement sequence of all sequences in a file and \n\t\t\toutputs them to a new file. Can handle all IUPAC ambiguity codes.\n\t\tUsage\n\t\t\tpython seqTools_v2.3.py batchRevComp foo.fasta [options]\n\t\tOptions\n\n\tbatchTranslate\n\t\t\tTool to translate many fasta files. Takes a file of file names as input and \n\t\t\tsequenctially opens each file and translates it by calling the buildCodonDict \n\t\t\tsubmodule.\n\t\tUsage\n\t\t\tpython seqTools_v2.3.py batchTranslate foo.fofn [options]\n\t\tOptions\n\n\tvcf2Fasta\n\t\t\tTool that produces a fasta sequence from a Variant Call Formatted (VCF) file. Tool uses \n\t\t\tstandard IUPAC ambiguity codes. Current filters include minimum \n\t\t\tdepth over the site and mapping quality assigned by the mapping program. \n\t\tUsage\n\t\t\tpython seqTools_v2.3.py vcf2Fasta foo.vcf [options]\n\t\tOptions\n\t\t\t-rc\tReverse complement. Output will be the reverse complement of the mapped sequence.\n\t\t\t\t[False]\n\t\t\t\tTrue, False\n\t\t\t-d\tMinimum depth requirement. Base calls with depth of coverage lower than -d will be \n\t\t\t\toutput as missing data (N). [10]\n\t\t\t-q\tMinimum mapping quality threshold. Base calls with lower Phred-scaled mapping \n\t\t\t\tquality than -q will be output as missing data (N). [30]\n\n\tseqRenamer\n\t\t\tTool to rename sequences within a fasta file. Input must be fasta. Output will be another \n\t\t\tfasta file with "renamed.fasta" as its suffix, and the original file name as its prefix.\n\t\t\tRenaming scheme can be numeric (default) or a fixed name. Numeric renaming will add "_n"\n\t\t\tto the end of the sequence tag, where n = the order of the sequence in the input fasta. \n\t\t\tFixed renaming will add a user-specified string to the beginning of each sequence tag \n\t\t\tlike follows: "><STRING>_sequenceTag".\n\t\tUsage\n\t\t\tpython seqTools_v2.3.py seqRenamer foo.fasta [options]\n\t\tOptions\n\t\t\t-n\tRenaming method. If the default numeric renaming system is not used, the value\n\t\t\t\tfollowing the option flag will be appended as a prefix to all sequences in the \n\t\t\t\tfasta file. [numeric]\n\n\tfasta2fastq\n\t\t\tMakes a fastq file with fake quality scores (ASCII_BASE 33) of -q from a fasta input. \n\t\t\tNote: input file handle MUST BE ".fasta"\n\t\tUsage\n\t\t\tpython seqTools_v2.3.py fasta2Fastq foo.fasta [options]\n\t\tOptions\n\t\t\t-q\tPhred quality score for fake quality scores. Must be an integer. [40]\n\n\tnexus2Fasta\n\t\t\tTool for converting nexus files into fasta files. Input should be a nexus file. Output \n\t\t\twill be a fasta file.\n\t\tUsage\n\t\t\tpython seqTools_v2.3.py nexus2Fasta foo.nexus [options]\n\t\tOptions\n\n\tfasta2Nexus\n\t\t\tTool for converting fasta files to nexus format. Input should be a fasta file, output \n\t\t\twill be a nexus file. Nexus file will be compatible with Paup or MrBayes, depending on \n\t\t\toption selected. Additional compatibility can be added by contacting the developer.\n\t\tUsage\n\t\t\tpython seqTools_v2.3.py fasta2Nexus foo.fasta [options]\n\t\tOptions\n\t\t\t-p\tIdentifies the compatibility mode desired. If you would like additional compatibility \n\t\t\t\tadded, please email developer. If paup compatibility is specified, a paup block at \n\t\t\t\tthe bottom of the nexus file will be written, which instructs paup to construct a \n\t\t\t\tNeighbor Joinging Tree and a ML pairwise distance matrix of the alignment. [none]\n\t\t\t\tnone,mrbayes, paup\n\n\thelp\n\t\t\tWrites this menu to the standard output.\n\t\tUsage\n\t\t\tpython seqTools_v2.3.py help \n\t\t\t\nInput should be in the form required by the submodule command. Note the program is case sensitive.\n\nIf you would like to know which command to use, type:\n\n\tpython seqTools_v2.3.py help\n\nQuestions should be directed to Joel Sharbrough at jsharbro[at]gmail.com\n\nCopyright (c) 2017 Joel Sharbrough'


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
        elif nuc == 'V':
            seq_revc = 'B' + seq_revc
        elif nuc == 'H':
            seq_revc = 'D' + seq_revc
        elif nuc == 'D':
            seq_revc = 'H' + seq_revc
        elif nuc == 'B':
            seq_revc = 'V' + seq_revc
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
    ambiguityCodes = {"A/C":"M","A/G":"R","A/T":"W","C/G":"S","C/T":"Y","G/T":"K","C/A":"M","G/A":"R","T/A":"W","G/C":"S","T/C":"Y","T/G":"K",'A/C/G':'V','C/A/G':'V','A/G/C':'V','C/G/A':'V','G/A/C':'V','G/C/A':'V','A/C/T':'H','C/A/T':'H','A/T/C':'H','C/T/A':'H','T/A/C':'H','T/C/A':'H','A/T/G':'D','T/A/G':'D','A/G/T':'D','T/G/A':'D','G/A/T':'D','G/T/A':'D','T/C/G':'B','C/T/G':'B','T/G/C':'B','C/G/T':'B','G/T/C':'B','G/C/T':'B',"A/T/C/G":"N","A/T/G/C":"N","A/C/T/G":"N","A/C/G/T":"N","A/G/C/T":"N","A/G/T/C":"N","T/A/C/G":"N","T/A/G/C":"N","T/C/A/G":"N","T/C/G/A":"N","T/G/C/A":"N","T/G/A/C":"N","C/A/T/G":"N","C/A/G/T":"N","C/T/A/G":"N","C/T/G/A":"N","C/G/A/T":"N","C/G/T/A":"N","G/A/C/T":"N","G/A/T/C":"N","G/C/A/T":"N","G/C/T/A":"N","G/T/A/C":"N","G/T/C/A":"N"}
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
                refSplit = ref.split(',')
                nucleotideList = ''
                for item in refSplit:
                    nucleotideList+= item + '/'
                nucleotideList = nucleotideList[0:-1]
                if '/' not in nucleotideList:
                    nucleotide = ref
                elif nucleotideList in ambiguityCodes:
                    nucleotide = ambiguityCodes[nucleotideList]
                else:
                    nucleotide = 'N'
            elif gt == '1/1':
                altSplit = alt.split(',')
                nucleotideList = ''
                for item in altSplit:
                    nucleotideList+= item + '/'
                nucleotideList = nucleotideList[0:-1]
                if '/' not in nucleotideList:
                    nucleotide = alt
                elif nucleotideList in ambiguityCodes:
                    nucleotide = ambiguityCodes[nucleotideList]
                else:
                    nucleotide = 'N'
            elif gt == '0/1':    
                refSplit = ref.split(',')
                altSplit = alt.split(',')
                nucleotideList = ''
                for item in refSplit:
                    nucleotideList += item + '/'
                for item in altSplit:
                    nucleotideList += item + '/'
                nucleotideList = nucleotideList[0:-1]
                if nucleotideList in ambiguityCodes:
                    nucleotide = ambiguityCodes[nucleotideList]
                else:
                    nucleotide = 'N'
            elif gt == '1/2':
                altSplit = alt.split(',')
                nucleotideList = ''
                for item in altSplit:
                    nucleotideList += item + '/'
                nucleotideList = nucleotideList[0:-1]
                if nucleotideList in ambiguityCodes:
                    nucleotide = ambiguityCodes[nucleotideList]
                else:
                    nucleotide = 'N'
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

def seqRenamer(fasta,renameMethod='numeric'):
    seqDict,seqList = buildSeqDict(fasta)
    outfile = open(fasta[0:-6] + '_renamed.fasta','w')
    if renameMethod == 'numeric':
        i = 1
        for seq in seqList:
            outfile.write(seq + '_' + str(i) + '\n' + seqDict[seq] + '\n')
            i += 1
    else:
        for seq in seqList:
            outfile.write('>' + renameMethod + '_' + seq[1:] + '\n' + seqDict[seq] + '\n')
    outfile.close()
    
def fasta2Fastq(fasta,qual='40'):
    #Makes a fastq file with fake quality scores of qual [default 40] from a fasta input. Note: input file handle MUST BE ".fasta"
    qualDict = {0:'!',1:'"',2:'#',3:'$',4:'%',5:'&',6:"'",7:'(',8:')',9:'*',10:'+',11:',',12:'-',13:'.',14:'/',15:'0',16:'1',17:'2',18:'3',19:'4',20:'5',21:'6',22:'7',23:'8',24:'9',25:':',26:';',27:'<',28:'=',29:'>',30:'?',31:'@',32:'A',33:'B',34:'C',35:'D',36:'E',37:'F',38:'G',39:'H',40:'I',41:'J',42:'K'}
    qualVal = qualDict[int(qual)]
    fileName = fasta[0:-5] + 'fastq'
    outfile = open(fileName,'w')
    seqDict,seqList = buildSeqDict(fasta)
    for gene in seqList:
        outfile.write('@' + gene[1:] + '\n')
        outfile.write(seqDict[gene] + '\n')
        outfile.write('+\n')
        outfile.write(qualVal*(len(seqDict[gene]))+'\n')
    outfile.close()

def nexus2Fasta(nexus):
    seqList = []
    seqDict = {}
    outfile = open(nexus[0:-3] + 'fasta','w')
    infile = open(nexus,'r')
    taxa = False
    characters = False
    seq = False
    currSeq = ''
    seqTag = ''
    for line in infile:
        newLine = line
        if len(newLine) > 3:
            while newLine[-1] == '\n' or newLine[-1] == '\t' or newLine[-1] == '\r':
                newLine = newLine[0:-1]
            if newLine[0:11] == 'BEGIN TAXA;':
                taxa = True
            elif taxa == True and newLine[0] == "'":
                lineSplit = line.split("'")
                seqTag = lineSplit[1]
                while seqTag[-1] == ';':
                    seqTag = seqTag[0:-1]
                seqList.append(seqTag)
            elif newLine[0:18] == 'BEGIN CHARACTERS;':
                characters = True
            elif characters == True and newLine[0] == "'":
                seq = True
                if currSeq != '':
                    seqDict[seqTag] = currSeq
                newLineSplit = newLine.split("'")
                seqTag = newLineSplit[1]
                currSeq = ''
                charNum = len(seqTag) + 2
                for char in newLine[len(seqTag):]:
                    if char == ' ':
                        charNum += 1
                currSeq = newLine[charNum:]
                while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r' or currSeq[-1] == ';':
                    currSeq = currSeq[0:-1]
            elif newLine == 'END;':
                taxa = False
                if characters == True and seq == True:
                    seqDict[seqTag] = currSeq
                characters = False
                seq = False
            elif characters == True and seq == True:
                currSeq += newLine
                while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r' or currSeq[-1] == ';':
                    currSeq = currSeq[0:-1]
    refSeqName = seqList[0]
    refSeq = seqDict[refSeqName]
    for seq in seqList:
        i = 0
        currSeq = seqDict[seq]
        outfile.write('>' + seq + '\n')
        for char in currSeq:
            if char == '.':
                outfile.write(refSeq[i])
            elif char == '?':
                outfile.write('N')
            else:
                outfile.write(char)
            i += 1
        outfile.write('\n')
    infile.close()
    outfile.close()

def fasta2nexus(fasta,program='none'):
    paupStatement = 'begin paup;\n\tnj treefile=N-mt_NJ.tre brlens=yes append=yes;\n\tlscore 1/nst=6 rmatrix=estimate rates=gamma shape=estimate;\n\tlset rmatrix=previous shape=previous;\n\tdset distance=ml missdist=ignore;\n\tsavedist NDecimals=4 file=' + fasta[0:-6] + '_dist.txt replace=yes undefined=asterisk;\nend;'
    outfile = open(fasta[0:-5] + 'nexus','w')
    seqDict,seqList = buildNexusSeqDict(fasta)
    longestSeqName = 0
    for seq in seqList:
        if len(seq) > longestSeqName:
            longestSeqName = len(seq)
    if program == 'mrbayes':
        outfile.write('#NEXUS\n\n\nbegin data;\n\tdimensions ntax=' + str(len(seqDict)) + ' nchar=' + str(len(seqDict[seqList[0]])) + ';\n\tformat datatype=dna missing=? gap=-;\n\tmatrix\n\t\t')
        for taxon in seqList:
            outfile.write(taxon + (' '*(longestSeqName - len(seq))) + '\t\t' + seqDict[taxon] + '\n\t\t')
        outfile.write(';\nend;\n')
    elif program == 'paup':
        outfile.write('#NEXUS\n\n\nbegin taxa;\n\tdimensions ntax=' + str(len(seqDict)) + ';\n\ttaxlabels')
        for seq in seqList:
            outfile.write(' ' + seq)
        outfile.write(';\nend;\n\nbegin characters;\n\tdimensions nchar=' + str(len(seqDict[seqList[0]])) + ';\n\tformat datatype=dna missing=? gap=- labels=left;\n\tmatrix\n')
        for taxon in seqList:
            outfile.write('\t\t' + taxon + (' '*(longestSeqName - len(seq))) + '\t\t' + seqDict[taxon] + '\n')
        outfile.write('\t\t;\nend;\n\n')
        outfile.write(paupStatement)
    else:
        outfile.write('#NEXUS\n\n\nbegin taxa;\n\tdimensions ntax=' + str(len(seqDict)) + ';\n\ttaxlabels ')
        for seq in seqList:
            outfile.write(' ' + seq + ';\nend;\n\nbegin characters;\n\tdimensions nchar=' + str(len(seqDict[seqList[0]])) + ';\n\tformat datatype=dna missing=? gap=- labels=left;\n\tmatrix\n\t\t')
        for taxon in seqList:
            outfile.write(taxon + (' '*(longestSeqName - len(seq))) + '\t\t' + seqDict[taxon] + '\n\t\t')
        outfile.write(';\nend;\n')
    outfile.close()

def buildNexusSeqDict(fasta):
    infile = open(fasta,'r')
    scaffoldDict = {}
    scaffoldList = []
    seqName = ''
    currSeq = ''
    for line in infile:
        if line[0] == '>':
            if seqName != '':
                newSeq = ''
                for char in currSeq:
                    if char == 'N' or char == 'n':
                        newSeq += '?'
                    else:
                        newSeq += char
                scaffoldDict[seqName] = newSeq
            seqName = line[1:]
            seqNameSplit = seqName.split('_')
            seqName = seqNameSplit[0]
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

def help():
    sys.stdout.write(readMe)

def buildCommand():
    command = sys.argv[1]
    programList = ['buildCodonDict','buildSeqDict','reverseComplement','cigarCalc','unwrap','batchRevComp','batchTranslate','vcf2Fasta','seqRenamer','help']
    programDict = {'batchRevComp': [], 'batchTranslate': [], 'buildCodonDict': ['-c'], 'buildSeqDict': [], 'cigarCalc': [], 'help': [], 'reverseComplement': [], 'unwrap': [], 'vcf2Fasta': ['-rc', '-d', '-q'],'seqRenamer':['-n'],'fasta2Fastq':['-q']}
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
                logfile.write(item + ' ' + str(optionValues[j]) + '\t')
                j += 1
            logfile.write('\n')
            logfile.close()
            commandStatement = command + '("' + Input + '","' + str(optionValues[0]) + '","' + str(optionValues[1]) + '","' + str(optionValues[2]) + '")'
        elif command == 'seqRenamer':
            optionValues = ['numeric']
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
        elif command == 'fasta2Fastq':
            optionValues = ['40']
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
        elif command == 'fasta2Nexus':
            i = 0
            optionValues = ['none']
            while i+4 < len(sys.argv):
                if sys.argv[i+3] == '-p':
                    if sys.argv[i+4] != 'none':
                        optionValues[0] = sys.argv[i+4]
                i += 1
            logfile = open('seqTools.log','a')
            logfile.write('Running ' + command + ' assuming the following parameters:\n\t')
            j = 0
            for item in options:
                logfile.write(item + ' ' + str(optionValues[j]) + '\t')
                j += 1
            logfile.write('\n')
            logfile.close()
            commandStatement = command + '("' + Input + '","' + str(optionValues[0]) + '","' + str(optionValues[1]) + '","' + str(optionValues[2]) + '")'
    exec(commandStatement)
    
buildCommand()

