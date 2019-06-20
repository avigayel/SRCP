# SRCP
## Short Read CircRNA Pipeline

### 1 Overview
The Short Reads circRNA Pipeline (SRCP) is made for identification of circRNAs when the length
of the reads is relatively small (say <75bp)
In such short reads the regular PCP pipeline identifies much less circle reads the SRCP.

### 2 Parameters
-circ (also –c or -C) – A file containing a list of circles locations in bed6 format.
-annotation (also –a or -A) – An annotation in REAL bed format
<chr><locA><locB><name><reads><strand><cdsStart><cdEnd><numExons><exon starts><exon ends>
-data (also –d or -dat) – the dataset file names separated by commas (explicitly in the shell)
Example: -d <file1>,<file2>,...<fileN>

### 3 Optional
-read_len (also –len or -l) – Explicitly give the script the size of the index it should create.
DEFAULT = automatically calculate for each sample based on reads length
-junction (also –j or -junc) – The length of the junction flanking region. DEFAULT = 10
-system (also –s or -sys) – The system genome to use. Choose from: dm3,mm9,hg19,rn6 or
the directory with the system with the prefix. i.e. /home/kadenerlab/Documents/reference/dm3/dm3
The last dm3 is not a directory but rather the prefix of all index and .fa files for the system, the same format that bowtie takes its indexes. Should contain also the .fa file with genome. DEFAULT = dm3
-transcriptome (also –trans or -tra) – The system transcriptome to use.
The directory with the system with the prefix. i.e. /home/kadenerlab/Documents/reference/riboIP_con_transcriptome/myGTFdm3_multiFasta_all
The last 'myGTFdm3_multiFasta_all' is not a directory but rather the prefix of all index files for the system, the same format that bowtie takes its indexes.
DEFAULT =/home/kadenerlab/Documents/reference/riboIP_con_transcriptome/myGTFdm3_multiFasta_all
-numThread (also -p) – The number of threads to use when possible. DEFAULT = 1
-pcp – Custom linear junctions, will be added to the index created for the linear junctions.
-mis (also -m or -M) – The maximal edit distance allowed for the index alignment in circular junctions and the maximal edit distance allowed in all alignments for linear junctions. DEFAULT=2
-N – Allow N mismatch in the bowtie alignment. Should only be used with reads of length < 32.
-safety (also -safe or -sa) – The minimal allowed safety margin, calculated as:
min(edit_distance_genome-edit_distance_index , edit_distance_transcriptome – edit_distance_index)
DEFAULT = maximal edit distance + 1
-quality (also -qual or -q) – The quality standard . DEFAULT = Phred+33
-min_quality (also -min_qual or -mq) – The minimum quality of the read (average) below which reads would be discarded. DEAFULT = 30
-output (also –o or -O) – The output file name. DEFAULT = “outputIndex”
--only_annotated – Search for only annotated(EXACT) circles. No internals. DEAFULT = false
--keep_temporary (also --keep) – Keep temporary runtime files. DEFAULT=false
(otherwise removes tmp/ directory)
--version (also –v) – Show the version and about.

### 4 Remarks
The final output file containing the: num unique circ reads, num circ reads,num linear reads left, num unique left linear reads, num right linear reads, num unique right linear reads ratios.
<output filename>_allSample.txt
The script creates a lot of other intermediate files. The circ and leftLinear and rightLinear circles with the junction are also output where the script was initiated. Other files are in the tmp/ directory.
