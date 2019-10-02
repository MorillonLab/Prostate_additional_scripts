# Prostate_additional_scripts
Additional scripts used for the prostate paper "Blind exploration of the unreferenced transcriptome reveals novel RNAs for prostate cancer diagnosis"

---

## 1) Differential expression (gencode27lift37 + holdUP annotation) ##

---


- 1-a) Design :

|condition| sample name |
|--------  |-------------|
|tumor	|	B67T18|
|tumor		|	B67T19|
|tumor		|	B67T20|
|tumor		|	B67T22|
|tumor		|	B65T1|
|tumor		|	B65T3|
|tumor		|	B65T5|
|tumor		|	B65T7|
|tumor		|	B67T11|
|tumor		|	B67T17|
|tumor		|	B67T21|
|tumor		|	B67T12|
|tumor		|	B67T23|
|tumor		|	B67T24|
|tumor		|	B67T9|
|tumor		|	B67T16|
|normal		|			B65T2|
|normal			|		B65T4|
|normal			|		B65T6|
|normal			|		B65T8|
|normal			|		B67T10|
|normal				|	B67T13|
|normal			|		B67T14|
|normal				|	B67T15|

```diff
- remark : the fastq files can be found in the GEO accession of the paper, as well as the instructions for the alignment !
```

 - 1-b) Get counts (gencode annoation + holdUP) : 
 
   - script (adapt the "input data" part) : https://github.com/MorillonLab/Prostate_additional_scripts/blob/master/prostate_counting.sh
   
   - link for gencode annotation (variable "gff_gtf_ForCounting" in the script) : ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.annotation.gff3.gz (uncompress it)
   
   ```diff
   - remark : the count files can also be found in the GEO accession of the paper !
   ```
   
   
 - 1-c) Comparison of conditions "normal" vs "tumoral" : 
 
   - script : https://github.com/MorillonLab/Prostate_additional_scripts/tree/master/DESeq2_prostate24
   
        - This directory contains the script `DESeq2_script.R` and its readme `README.txt`
        
   - How to use it : read the readme `README.txt`, or in a terminal, type : ./DESeq2_script.R to see the usage
   
   ---
   
 ## 2) Overlap between GENCODE, MiTranscriptome, DE-kupl and HoLdUp catalogues ##
   
   ---
   
  - 2-a) script (adapt the "input data" part) : https://github.com/MorillonLab/Prostate_additional_scripts/blob/master/getRecursiveIntersections.sh
  
      -> check the results in the table `all_intersections.tsv`
   
   ---
   
 ## 3) Unsupervised clustering of prostate specimens
     
   ---
   
   - 3-a) script (adapt the "input data" part) : https://github.com/MorillonLab/Prostate_additional_scripts/blob/master/dekupl_heatmap_venn.R
   
      -> replace the gff files by :
      
      - https://github.com/MorillonLab/Prostate_additional_scripts/blob/master/data/kmers_contiguous.gff
      - https://github.com/MorillonLab/Prostate_additional_scripts/blob/master/data/kmers_repeat.gff
      - https://github.com/MorillonLab/Prostate_additional_scripts/blob/master/data/kmers_spliced.gff
      - https://github.com/MorillonLab/Prostate_additional_scripts/blob/master/data/kmers_unmapped.gff
           
