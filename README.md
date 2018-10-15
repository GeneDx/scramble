Scramble
========

Repository Layout
-----------------
1. cluster_identifier/ - this directory contains the part of the application responsible for identifying soft clipped clusters. For how to build see the build section.
2. cluster_analysis/ - contains code related to the analysis of soft-clipped clusters and interpretation.
3. validation/ - sample bam and outputs for testing SCRAMble
4. README - Repository description/notes.

Build
-----
To build the cluster_identifier:
$ cd cluster_identifier/src
$ make

That should be it. It will create an executable named cluster_identifier. 


Running
-------
SCRAMBle runs as a two-step process. First cluster_identifier is used to generate soft-clipped read cluster consensus sequences.
Second, SCRAMBle-MEIs.R analyzes the cluster file for likely MEIs.

To run SCRAMble cluster_identifier:
/path/to/install_dir/scramble/cluster_identifier/src/build/cluster_identifier /path/to/install_dir/scramble/validation/test.bam

To run SCRAMBle-MEIs (with default settings):
Rscript --vanilla /path/to/install_dir/scramble/cluster_analysis/bin/SCRAMble-MEIs.R \
	-o /path/to/output/out.txt 	\
	-c /path/to/install_dir/scramble/validation/test.clusterFile.txt \
	-i /path/to/install_dir/scramble/cluster_analysis/bin/ \
	--mei-refs /path/to/install_dir/scramble/cluster_analysis/resources/MEI_consensus_seqs.fa  
	

Output
------
The output of cluster_identifier is a tab delimited text file with clipped cluster consensus sequences.
The columns are as follows:
	1. Coordinate
	2. Side of read where soft-clipped ocurred
	3. Number of reads in cluster
	4. Clipped read consensus
	5. Anchored read consensus
	
The output of SCRAMBle-MEIs.R is a tab delimited text file with MEI calls. If no MEIs are present an output file will still be produced with only the header.
The columns are as follows:
	1. Insertion 					=> coordinate where MEI insertion occurs (zero-based)
	2. MEI_Family 					=> the consensus sequence to which the clipped sequence aligned best
	3. Insertion_Direction 			=> whether MEI is on fwd or rev strand relative to bam reference
	4. Clipped_Reads_In_Cluster		=> number of supporting reads in cluster
	5. Alignment_Score 				=> pairwise alignment score of clipped read consensus to MEI reference sequence
	6. Alignment_Percent_Length 	=> percent of clipped read consensus sequence involved in alignment to MEI reference sequence
	7. Alignment_Percent_Identity 	=> percent identify of alignment of clipped read consensus sequence with MEI reference sequence
	8. Clipped_Sequence				=> clipped cluster consensus sequences
	9. Clipped_Side					=> left or right, side of read where soft-clipping ocurred
	10. Start_In_MEI				=> left-most position of alignment to MEI reference sequence
	11. Stop_In_MEI					=> right-most position of alignment to MEI reference sequence
	12. polyA_Position				=> position of  polyA clipped read cluster if found
	13. polyA_Seq					=> clipped cluster consensus sequences of polyA clipped read cluster if found
	14. polyA_SupportingReads		=> number of supporting reads in polyA clipped read cluster if found
	15. TSD							=> target site duplication sequence if polyA clipped read cluster found
	16. TSD_length					=> length of target site duplication if polyA clipped read cluster found

	
	
R Dependencies
--------------
SCRAMBLE-MEIs.R was developed on R version 3.1.1 (2014-07-10) and uses the following libraries:
[1] Biostrings_2.34.1   XVector_0.6.0       IRanges_2.0.1
[4] S4Vectors_0.4.0     BiocGenerics_0.12.1 stringr_0.6.2
[7] optparse_1.4.4


Disclaimers
-----------
In theory, SCRAMble should work well on any MEI reference fasta sequences, however, it has only been tested on the sequences provided in /path/to/install_dir/scramble/cluster_analysis/resources/MEI_consensus_seqs.fa


