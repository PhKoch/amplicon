# Description
The script `identifyAmplicons.py` searches sequences of amplicons stored in a fastq file based on primer sequences that allow to identify the amplicon.

# Usage

## Command

`identifyAmplicons.py [-h] -p PRIMERSFILE [-d] fastq`

* `identifyAmplicons.py -h` for details on the command line options.

* `identifyAmplicons.py --primers primers.fa Pool_A.merged.800.fastq` to find matches of some test primers in a fastq file of 200 merged reads.

## Input

### `PRIMERSFILE`
A fasta file storing the pairs of primers in the following format:

```
>Amplicon1Fwd Amplicon1
GGAAGAG
>Amplicon1Rev Amplicon1
CCAGCGC
>Amplicon2Fwd Amplicon2
...
>Amplicon2Rev Amplicon2
...
```

* First entry: forward primer in 3' - 5' direction
* Second entry: reverse primer in 3' - 5' direction
* Amplicons will be paired, based on the second field of the identifier line. Therefore, multiple pairs in one file are allowed.
* See example in `primers.fa`


### `fastq`

A normal fastq file, for example produced on an Illumina HiSeq or MiSeq machine. **Adapters need to be removed** from the ends of the reads.

# Algorithm
The "Algorithm" is pretty basic and can be described as below.
```
-for each read
	-for each primer pair
		-look wether the read starts and ends with the repective primer sequence
		-look wether the read starts and ends with the reverse complement of the repective primer sequence
		-store match and read length in a file
-store unmatched reads in a file
```

# Notes & Restrictions

* The script was tested on and applied to 2x300bp paired-end MiSeq reads that were merged with `usearch`, as some expected amplicons were `>300bp`. In principle, it should work also on single files (R1/R2, one per run) of a paired end data set.
* By default, amplicons are searched in both directions, so reverse complements will automatically be calculated.
* If the reverse complement of a primer pair identifies an amplicon on the negative strand, its output in the `matched.csv` file will be reverse complemented to have unified sequences. Ouf course, such a hit on the minus strand is marked with a `-` in the `read_strand` field.


# Compatibility
tested on Python 3.7.2 & 3.8.2
