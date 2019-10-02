#!/bin/bash


#### input data ########

#path to bedtools program
bedtools="bedtools"

#path to samtools program
samtools="samtools"


#list of annotations
#the annoation of the contigs (kmers_contiguous.gff, kmers_spliced.gff...) can be found in the directory data (https://github.com/MorillonLab/Prostate_additional_scripts/tree/master/data)
#gencode gene level annotation can be found here : https://github.com/MorillonLab/Prostate_additional_scripts/blob/master/data/gencode.v27lift37_gene_lvl.gff
#MiTranscriptome lncRNAs can be found here (uncompress it) : https://github.com/MorillonLab/Prostate_additional_scripts/blob/master/data/mitranscriptome.v2_lncRNAs_gffread.zip
annot_list=("/home/marcgabriel/Documents/gencode27lift37/gencode.v27lift37_gene_lvl.gff"
            "/home/marcgabriel/Documents/Marina_P/mitranscriptome.gtf/mitranscriptome.v2_lncRNAs_gffread.gff"
            "/home/marcgabriel/Documents/Marina_P/data_intersections/holdup_prostate_class1w_only.gff"
	    "/home/marcgabriel/Documents/Marina_P/kmers_annotations/kmers_contiguous.gff"
            "/home/marcgabriel/Documents/Marina_P/kmers_annotations/kmers_spliced.gff"
            "/home/marcgabriel/Documents/Marina_P/kmers_annotations/kmers_repeat.gff"
            "/home/marcgabriel/Documents/Marina_P/kmers_annotations/kmers_unmapped.gff")
            

#names associated to the annotations (order respected)
prefix=("gencode"
        "MiTranscriptome_lncRNAs"
        "holdUP_all_class1"
        "dekupl_contigs_contiguous"
        "dekupl_contigs_spliced"
        "dekupl_contigs_repeat"
        "dekupl_contigs_unmapped")

#output directory            
output_dir="/home/marcgabriel/Desktop/test4"


####################

output_dir="${output_dir}/"
output_dir=$(echo $output_dir |sed 's/\/\//\//g')
if [ ! -d $output_dir ]; then mkdir -p $output_dir; fi
            

#put the header in the final table
echo -e "cond1\tcond2\tcond1_counts\tcond2_counts\tminimal_overlap\tmax_overlap_cond1\tmax_overlap_cond2\tmax_overlap_cond1_percent\tmax_overlap_cond2_percent" >${output_dir}all_intersections.tsv


########## intersection with other databases ##############

#percentage of reciprocal overlap (bedtools parameter)
#overlap_frac=""
overlap_frac="-f 0.50 -F 0.50 -e"
#overlap_frac="-f 0.20 -F 0.20 -e"


echo -e "fraction of overlap between annotation 1 & 2 : $overlap_frac\n" >${output_dir}summary.txt

summary="${output_dir}summary.txt"

#all possible combinations in the annot list (for 3 annotations for example, we have 1 vs 2, 1vs 3, & 2 vs 3)
#in bash, incrementation in lists start from 0, contrary to R
#pairwise intersections

#index example : 0,1,2,3...
for i in $(seq 0 $((${#annot_list[*]}-1)));do

  #shift the index table
  #new index example : 1,2,3,4...
  new_list=($(seq $((i+1)) $((${#annot_list[*]}-1))))


  for j in ${new_list[*]};do
  
  
         
       
	  cond1_counts=$(wc -l ${annot_list[$i]} |awk '{print $1}')
	  
	  cond2_counts=$(wc -l ${annot_list[j]} |awk '{print $1}')
	  
	  echo -e "- annotations to intersect : ${annot_list[$i]} & ${annot_list[j]} (= ${prefix[$i]} & ${prefix[j]})\n" >>$summary
	  
	  echo -e "type1 : $type1; type2 : $type2\n" >>$summary
	  
	  echo -e "\t- number of feat in ${prefix[$i]} : $cond1_counts\n" >>$summary
	  
	  echo -e "\t- number of feat in ${prefix[j]} : $cond2_counts\n" >>$summary

	  $bedtools intersect $overlap_frac -wb -s -nobuf -s -nonamecheck -a ${annot_list[$i]} -b ${annot_list[j]} >${output_dir}intersected.tmp
	 
	  max_overlap_cond1=$(cut -f9 ${output_dir}intersected.tmp |grep -v "^$"|sort -u |wc -l )
	    
	    
	  max_overlap_cond2=$(cut -f18 ${output_dir}intersected.tmp |grep -v "^$"|sort -u|wc -l )
	  
	  
	  if [[ $max_overlap_cond1 -gt $max_overlap_cond2 ]];then
	  
	      shared_counts=$max_overlap_cond2
	      
	  else
	  
	      shared_counts=$max_overlap_cond1
	  
	  
	  fi
	  
	  percent_counts=($(awk -v cond1_counts=$cond1_counts -v cond2_counts=$cond2_counts -v shared_counts=$shared_counts 'BEGIN{OFS="\t";print (shared_counts/cond1_counts)*100,(shared_counts/cond2_counts)*100}'))
	  
	  percent_counts_diff=($(awk -v cond1_counts=$cond1_counts -v cond2_counts=$cond2_counts -v max_overlap_cond1=$max_overlap_cond1 -v max_overlap_cond2=$max_overlap_cond2 'BEGIN{OFS="\t";print (max_overlap_cond1/cond1_counts)*100,(max_overlap_cond2/cond2_counts)*100}'))
	  
	  percent_in_cond1=${percent_counts[0]}
	  
	  percent_in_cond2=${percent_counts[1]}
	  
	  percent_in_cond1_diff=${percent_counts_diff[0]}
	  
	  percent_in_cond2_diff=${percent_counts_diff[1]}
	  
	  echo -e "\t\t- number of common features : $shared_counts (= ${percent_in_cond1}% of ${prefix[$i]} & ${percent_in_cond2}% of ${prefix[j]})\n" >>$summary
	  
	  echo -e "\t\t- number of features for each (diff) : $max_overlap_cond1 = ${percent_in_cond1_diff}% of ${prefix[$i]} ; $max_overlap_cond2 = ${percent_in_cond2_diff}% of ${prefix[j]})\n" >>$summary
	  
	  echo -e "\n************\n" >>$summary
	  
	  
	  echo -e "${prefix[$i]}\t${prefix[j]}\t${cond1_counts}\t${cond2_counts}\t$shared_counts\t$max_overlap_cond1\t$max_overlap_cond2\t${percent_in_cond1_diff}\t${percent_in_cond2_diff}" >${output_dir}${prefix[$i]}_${prefix[j]}_intersections.tsv
	  
	  cat ${output_dir}${prefix[$i]}_${prefix[j]}_intersections.tsv >>${output_dir}all_intersections.tsv
	 
  
  done

done

echo -e "The results are in the directoty : ${output_dir}\n"
