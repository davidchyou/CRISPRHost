CRISPRHost.pl

GENERAL USAGE

	perl CRISPRHost.pl -in <fasta_path> -out <dir_path> <options>
	perl /path/CRISPRHost.pl -in <fasta_path> -out <dir_path> <options>

OPTIONS

	-in <path>                   [Required] A FASTA file of viral genome, mobile genetic elements, or sequences.
	-out <path>                  [Required] Path to the output directory.
	-kingdom <A|B|AB>            Select a kingdom (A=archaea, B=bacteria, AB=archaea+bacteria). 
	                             The corresponding BLASTDB of CRISPR spacers (non-redundant) will then be used:
  	                                1. A: RefSeq95 archaeal, DB/CRISPRBankSpacers_4_95_2555_100_archaea_refseq_nr.fa
    	                        2. B: RefSeq95 bacterial, DB/CRISPRBankSpacers_4_95_2555_100_bacteria_refseq_nr.fa
      	                        3. AB: RedSeq95 archaeal+bacterial, DB/CRISPRBankSpacers_4_95_2555_100_all_refseq_nr.fa
  	                             During the first run, the app will unzip "DB.zip" for the DB directory.
  	                             Default: AB.
	-mask_arrays                 Optionally run MINCED to predict predict arrays and mask them using BEDTOOLS.
	-dbsize <integer>            Number of nucleotides in BLASTDB, optional parameter for BLASTN. Default: 10000.
	-e_cutoff <numerical>        E-value cutoff for BLASTN spacer hits. Default: 0.001.
	-b_cutoff <numerical>        Bit-score cutoff for BLASTN spacer hits. Default: 50.
	-pr_cutoff <numerical 0-1>   Minimum proportion of base-paired nucleotides in spacer-target hybrids. Default: 0.7.
	-min_nocsh <integer>         Minimum number of CRISPR spacer hits per sequence. Default: 1.
	-max_target_seqs <integer>   Maximum number of hits to be returned by BLASTN for post-processing. Default: 500.
	-h                           Show this help.
	-help                        Show this help (same as -h).      

OUTPUT FILES

	summary.csv                Files showing the key results. The columns are:

                                   1. Spacer hit statistics
                             
                                   VIRAL_GENOME                Viral genome identifier.
                                   N_HOST_SPACER               Number of hits from spacers of the same sequence.
                                   E_VALUE                     BLASTN e-value.
                                   BITSCORE                    BLASTN bitscore.
                                   PROP_NT_MATCHED_IN_HYBRID   Proportion of matched nucleotide in a spacer-target hybrid.
                                
                                   2. Attributes of matched spacers
                             
                                   SPACER_ID                   Spacer identifier.
                                   KNOWN_SPECIES               A list of species where the same spacer can be found. 
                                                               These species are the predicted hosts.
                                   TARGET_ALIGN_START          Start-coordinate of the target alignment.
                                   SEQ_TARGET_ALIGNED          Sequence of the spacer target.
                                   SPACER_ORI_ALIGN            Whether the match is corresponding to the reverse-complement 
                                                               of the spacer (minus) or not (plus). 
                                                     
	full_results.csv             Files showing the key results in summary.csv, with additional information and statistics
                                 included.
                             
	spacer_matched.fna           A list of spacers matched in multi-FASTA format.

	c_arrays.gff                 Any CRISPR arrays found by MINCED and masked by BEDTOOLS, if -mask_arrays is used                         

DEPENDENCIES

1. NCBI BLAST suite: blastn, makeblastdb and blastdbcmd (provided)
2. Bedtools (provided)
3. minced (provided)
4. Perl library BioPerl
5. R 3.6+ (no add-on libraries needed) and Rscript

EXAMPLES

	perl CRISPRHost.pl -in NC_034623.fna -out test_NC_034623
	perl CRISPRHost.pl -in NC_034623.fna -out test_NC_034623 -kingdom AB
	perl CRISPRHost.pl -in NC_034623.fna -out test_NC_034623 -kingdom AB -mask_arrays
