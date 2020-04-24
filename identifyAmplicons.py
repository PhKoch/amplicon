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
	# "positional argument" for the fastq
	parser.add_argument('fastq', type=argparse.FileType('r', encoding='UTF-8'),
						help='Reads containing the amplicons in fastq format (required).')

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

	## read the primers file and get back structured records
	primerPairs = getPrimerPairs(args.primersfile)
	#print (primerPairs)

	# just a helper:
	linenumber = 1

	# prepare an output file
	# - the read sequence which is written to the end of a line is reverse complemented if the hit is on the minus strand
	# - by default we write only hits in this file. Can be extended later to write non hitting reads as well.
	outtabfile = open("outfile.txt",'w')
	# add the headerline to the file
	outtabfileHeaderline = 'read_id\t'
	outtabfileHeaderline += 'read_strand\t'
	outtabfileHeaderline += 'read_length\t'
	outtabfileHeaderline += 'amplicon_name\t'
	outtabfileHeaderline += 'left_primer_name\t'
	outtabfileHeaderline += 'left_primer_seq\t'
	outtabfileHeaderline += 'right_primer_name\t'
	outtabfileHeaderline += 'right_primer_seq\t'
	outtabfileHeaderline += 'read_seq\n'
	outtabfile.write(outtabfileHeaderline)

	# prepare an output file for not matching reads
	unmatchedreads = list()
	readnumber = 0
	matchnumber = 0
	# ATTENTION: some reads start with one "N" but would have hits. why is this???
	## iterate through the fqrecords
	for read in SeqIO.parse(args.fastq, "fastq"):
		readnumber = readnumber + 1
		#print ("%i %s" % (linenumber, read.id))
		## for each read we check all amplicons for a match,
		# so loop through the primer pairs
		foundmatch = False
		for amplicon in primerPairs:
			leftseq = primerPairs[amplicon][0]
			rightseq = primerPairs[amplicon][1]
			# check for match on + strand
			if(read.seq.startswith(leftseq.seq) & read.seq.endswith(rightseq.seq)):
				#print ("%i +++ %s" % (linenumber, read.id))
				#print (str(linenumber) + '     match for: ' + read.seq)
				#print (str(linenumber) + "     read length: %i | primer left: %s (%s) | primer right %s (%s)\n" % (len(read.seq), leftseq.id,leftseq.seq,rightseq.id,rightseq.seq))
				foundmatch = True
				matchnumber = matchnumber + 1
				outtabfile.write("%s\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.id,
																			'+',
																			len(read.seq),
																			amplicon,
																			leftseq.id,
																			leftseq.seq,
																			rightseq.id,
																			rightseq.seq,
																			read.seq))
			# check for match on - strand
			elif(read.seq.startswith(rightseq.seq.reverse_complement()) & read.seq.endswith(leftseq.seq.reverse_complement())):
				#print ("%i --- %s" % (linenumber, read.id))
				#print (str(linenumber) + '     match for: ' + read.seq)
				#print (str(linenumber) + "     read length: %i | primer left: %s (revcompl) (%s) | primer right %s (revcompl) (%s)\n" % (len(read.seq), rightseq.id,rightseq.seq.reverse_complement(),leftseq.id,leftseq.seq.reverse_complement()))
				foundmatch = True
				matchnumber = matchnumber + 1
				outtabfile.write("%s\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.id,
																			'-',
																			len(read.seq),
																			amplicon,
																			leftseq.id,
																			leftseq.seq,
																			rightseq.id,
																			rightseq.seq,
																			read.seq.reverse_complement()))
			# no hit for this amplicon
			#else:
				#print (str(linenumber) + ' NO match for: ' + read.seq + "\n")

		linenumber = linenumber + 1
		if foundmatch == False:
			unmatchedreads.append(read)

	# close the file handle as it was opend automatically by ArgumentParser
	args.fastq.close()
	outtabfile.close

	# write the unmatched reads to a file
	SeqIO.write(unmatchedreads, "unmatched_reads.fq", "fastq")
	print ("%i reads processed, found %i matches, wrote %i unmatched reads to the file unmatched_reads.fq" % (readnumber,matchnumber,len(unmatchedreads)))
if __name__ == "__main__":
	main()
