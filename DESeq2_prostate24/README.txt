1) DESCRIPTION :

  The script "DESeq2_script.R" does by default pairwise comparisons of conditions in order to find differentially expressed genes using the DESeq2 package.
 
  Go to the last part 4) to see how to run the script on the supplied data (prostate data, 8 normal vs 16 tumoral).


2) REQUIRED DEPENDANCIES :

  - Unix-like terminal (debian, ubuntu...) with grep installed.

  - R version >= 3.4.3, with :

		-a) R packages from CRAN
		
				ggplot2
				RColorBrewer
				optparse
				gplots
				pROC
				foreach
				doParallel	
				grid	
				BiocParallel
				data.table
		
		 -b) R packages from Bioconductor
	  
				 DESeq2
				 pheatmap

		 
3) USAGE :

   Type in the terminal (in the folder of the script) :
   
		  chmod 755 DESeq2_script.R 
		
   then :
   
	    ./DESeq2_script.R -f FILES_DESCRIPTOR -a GFF -o OUTPUT_DIR
	    

  "-f", "-a" and "-o" are the required parameters, the others are optional.
	
  Description of these parameters :
  
	  FILES_DESCRIPTOR : a tab-delimited file with 2 columns (full path to file + condition), like this :

		/path/to/file_1	condition_1
					  
		/path/to/file_2	condition_1
					  
		/path/to/file_3	condition_2
					  
		/path/to/file_4	condition_2
					  

	  GFF : annotation in gff3 format (uncompressed)

	  OUTPUT_DIR : output directory in which the results will be put.

  The results will be in "OUTPUT_DIR", and according to "FILES_DESCRIPTOR", a folder "condition_1_vs_condition_2" will be created.
  To see all the available options of the script, just run the script without any options (./DESeq2_script.R).


4) HOW TO RUN THE SCRIPT ON THE SUPPLIED DATA (in the sub-folder test_files)

  ./DESeq2_script.R -f ./test_files/counts_design.tsv -a ./test_files/gencode.v27lift37_sorted.gff3 -o test_results -n ./test_files/matching_gene_id_types.tsv

  In this example, we have 2 conditions (normal & tumoral, respectively with 8 & 16 samples) that will be compared.
  The results of this comparison will be in the folder "test_results".
  In the sub-folder "normal_vs_tumoral", you will find :
    
     - DESeq_output_normal_vs_tumoral.tsv           : raw output from DESeq2 (annotated with the gff file)
     - sig_diff_normal_vs_tumoral.tsv               : all significant differentially expressed genes according to the filters (default : padj <= 0.05)
     - sig_diff_downregulated_normal_vs_tumoral.tsv : significant downregulated genes (default : padj <= 0.05 & log2FC < 0)
     - sig_diff_upregulated_normal_vs_tumoral.tsv   : significant upregulated genes (default : padj <= 0.05 & log2FC > 0)
     - MAplots & volcano-plots for controls
     
     