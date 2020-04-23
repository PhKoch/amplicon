from Bio.Seq import Seq
from Bio import SeqIO

import argparse
import textwrap

import logging, sys

def pars():

	epilogText = """
	This is how PRIMERSFILE has to look like:

	>Amplicon1Fwd Amplicon1
	GGAAGAG
	>Amplicon1Rev Amplicon1
	CCAGCGC
	>Amplicon2Fwd Amplicon2
	...
	>Amplicon2Rev Amplicon2
	...

	First entry - forward primer.
	Second entry - reverse primer.
	Amplicons will be matched, based on the second field of the identifier line.
	Multiple paires are possible in the file.

	By default amplicons are searched in both directions, so reverse complements
	will automatically be calculated.
	"""

	parser = argparse.ArgumentParser(description='Finding amplicon sequences in a fastq file based on a list of primer sequences.',
									epilog = textwrap.dedent(epilogText),
									formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-p', '--primers', dest='primersfile', required=True,
						type=argparse.FileType('r', encoding='UTF-8'),
						help='Primer Sequences in fasta format (required). Details see below.')
	parser.add_argument('-d', '--debug', dest='loglevel', action='store_const',
						const=logging.DEBUG, default=logging.INFO, help='Show debug messages')

	args = parser.parse_args()

	# set the logging level (= having debug messages or not)
	# predefined levels and their integer values
	# CRITICAL 50
	# ERROR 40
	# WARNING 30
	# INFO 20
	# DEBUG 10
	# NOTSET 0
	logging.basicConfig(stream=sys.stderr, level=args.loglevel)

	#for debugging print the help:
	# parser.print_help()

	return parser, args


def getPrimerPairs(file):
	# read in the file of primers that we use to search the amplicons
	# in all input fastq files.
	primers = list(SeqIO.parse(file, "fasta"))
	# close the file handle as it was opend automatically by ArgumentParser
	file.close()

	primerPairs = dict()
	for record in primers:
		#print("%s %s %i" % (record.id, repr(record.seq), len(record)))
		ampName = record.description.split()[1]
		logging.debug("parsing entry for amplicon " + ampName + " " + record.seq)
		# if there is no entry for the ampName, build one
		if not primerPairs.get(ampName):
			primerPairs.update({ampName : record})
			logging.debug("create entry for " + ampName + " with sequence " + record.seq)
		# if there is already an entry for ampName, add the second sequences
		# (currently, this is not very robust as we do not ask for amplicons with >2 entries in the fasta file.)
		else:
			primerTupel = (primerPairs.get(ampName),record)
			primerPairs.update({ampName : primerTupel})
			logging.debug("added " + record.seq + " to " + ampName)
		#print (primerPairs)
		logging.debug(primerPairs[ampName])

	# print the full dictionary
	logging.debug(primerPairs)
	#print(primerPairs)

	# old output of the function (only the list of SeqRecord objects)
	# return primers

	# return a dictionary of structure {Amplicon_name : (SeqRecord,SeqRecord)}
	# that can be given to Seq.startswith()/endswith()
	return primerPairs


### here comes the main ###
def main():
	(parser, args) = pars()

	# #create a sequence object
	oneread = Seq("CCACACCCAGCTACTGACCTATACAACAGAAGGAAAAGTGCTATGTTTGGACATGGAATTCTTAGGAGCAACAGCTGGAAGAG")
	# #print out some details about it
	# print("seq %s is %i bases long" % (oneread, len(oneread)))
	# print("reverse complement is %s" % oneread.reverse_complement())
	# print("protein translation is %s" % oneread.translate())
	# print (dir(oneread))
	# if(oneread.startswith("CCACACCC") & oneread.endswith("GGAAGAG")):
	#	 print ('match for:' + oneread)

	## read the primers file and get back structured records
	primerPairs = getPrimerPairs(args.primersfile)
	#print (primerPairs)

	# test
	# for amplicon in primerPairs:
	#	 print(primerPairs[amplicon][1].id, primerPairs[amplicon][1].seq + ", revcomp: " + primerPairs[amplicon][1].seq.reverse_complement())


	# just a helper:
	linenumber = 1

	# read in the fastq file
	# Although it is elegant to use the iterator to go through the fastq
	# for read in SeqIO.parse("Pool_A.merged.800.fastq", "fastq"):
	# we turn the content into a list of SeqRecord items. Because we will
	# iterate multiple times through it and also want to calculate
	# the number of entries in the fastq file.
	fqrecords = list (SeqIO.parse("Pool_A.merged.800.fastq", "fastq"))

	# ATTENTION: some reads start with one "N" but would have hits. why is this???

	## loop through the fqrecords
	for read in fqrecords:
		#print ("%i %s" % (linenumber, read.id))
		## for each read we check all amplicons for a match,
		# so loop through the primer pairs
		for amplicon in primerPairs:
			# check for match on + strand
			if(read.seq.startswith(primerPairs[amplicon][0].seq) & read.seq.endswith(primerPairs[amplicon][1].seq)):
				print ("%i +++ %s" % (linenumber, read.id))
				print (str(linenumber) + '     match for: ' + read.seq)
				print (str(linenumber) + "     read length: %i | primer left: %s (%s) | primer right %s (%s)\n" % (len(read.seq), primerPairs[amplicon][0].id,primerPairs[amplicon][0].seq,primerPairs[amplicon][1].id,primerPairs[amplicon][1].seq))
			# check for match on - strand
			elif(read.seq.startswith(primerPairs[amplicon][1].seq.reverse_complement()) & read.seq.endswith(primerPairs[amplicon][0].seq.reverse_complement())):
				print ("%i --- %s" % (linenumber, read.id))
				print (str(linenumber) + '     match for: ' + read.seq)
				print (str(linenumber) + "     read length: %i | primer left: %s (revcompl) (%s) | primer right %s (revcompl) (%s)\n" % (len(read.seq), primerPairs[amplicon][1].id,primerPairs[amplicon][1].seq.reverse_complement(),primerPairs[amplicon][0].id,primerPairs[amplicon][0].seq.reverse_complement()))

			# no hit for this amplicon
			else:
				print (str(linenumber) + ' NO match for: ' + read.seq + "\n")
		linenumber = linenumber + 1
	# print (len(fqrecords))


if __name__ == "__main__":
	main()
