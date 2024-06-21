Run index.sh with all other scripts in the same directory.

Requires an indexed refrence genome (-R).  
Requires a directory with indexed .bam files (-B).  
Requires a .bed file listing genes of intrest (-G).  
Requires a list of .bam files and their sample names (-S).  
Optional .bed file listing exons of genes of intrest (-E).  
Optional number of bases to add to begining and end of gene (-p). Defaults to 10  
Optional depth of reads for a base to extract (-d). Defaults to 1  
Optional threshold for gene coverage (-c). Defaults to 90%  
Optional ID of previous job to contine from a previous job (-i)  
Optional step to continue from (-r). Allowed values: 'variant', 'extract', 'align', 'tree'.  

First step: Coverage  
	Submits an array job for each .bam listed in -S.  
 	Each job determines the coverage of all genes from -G (or just their exons if -E is provided).  
	It will attempt with 10 GB of RAM and if that fails, it will resubmit with addition resources up to 4 times (50 GB, 100 GB, 250 GB, 500 GB).  

index.sh will then determine which .bam files to use for each gene. By default, only .bams with over 90% coverage of a gene will be allowed to continue.  

Second step: Variant  
	Submits an array job for each gene listed in -G.  
 	Each job calls variant bases for each gene in all .bam files with enough coverage.  
	If the -S file contains 3 colums, the second column will be ued to determine groups for variant calling.  

 Third step: Extract  
 	Submits an array job for each gene listed in -G.  
	Each job extracts the gene of intrest from all .bam files with enough coverage.  
 	If depth of a read is to low, it will be read as 'N'.  
	For variants, IUPAC codes will be used.  

Fourth step: Align  
	Submits an array job for each gene listed in -G.  
	Each jobs aligns the various extracted genes with a copy of the gene from the refrence sequence.  
	It will attempt with 10 GB of RAM and if that fails, it will resubmit with addition resources up to 6 times (20 GB, 40 GB, 80 GB, 160 GB, 320 GB, 640 GB).  

 Fifth step: Tree
 	Submits an array job for each gene listed in -G.
  	Each job makes determining the best model to use and creates a phylogenetic tree of the previous alignment.
   	It will attempt with 20 GB of RAM and if that fails, it will resubmit with addition resources up to 5 times (40 GB, 80 GB, 160 GB, 320 GB, 640 GB).

At each step, the logs of succesful runs will be copied to the output directory and a tally of stats will be generated for each array.  

Finally, all the logs from the run will be moved to a single directory.  
