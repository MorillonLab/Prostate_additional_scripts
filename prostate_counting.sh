#!/bin/bash


#annotation
gff_gtf_ForCounting="/home/marcgabriel/Documents/gencode27lift37/gencode.v27lift37_sorted.gff3"

#it works like hseq-count : 1 -> read 1 is in the same strand as the RNA ; 2 -> read 1 is antisense to the RNA
orientation="-s 2"

#correspondance gencode gene types & global gene types
gene_type="/home/marcgabriel/Documents/gencode26/gene_types.tsv"

#featureCounts program
featureCounts_prog="/home/marcgabriel/Desktop/subread-1.6.0-source_MAX_HIT_NUMBER_10e6/bin/featureCounts"


#bedtools 
bedtools="bedtools"

#home-made annotation (from holdUp)

#class 1
init_holdupProstate1_class1w="/home/marcgabriel/Documents/Marina_P/lncRNA_gtf_HoLDuP_run_10/contigs_sans_sens.anti.class1.quantile.0.2.gtf"

init_holdupProstate2_class1w="/home/marcgabriel/Documents/Marina_P/lncRNA_gtf_HoLDuP_run_10/contigs_sans_sens.inter.class1.quantile.0.2.gtf"

#class 3
init_holdupProstate1_class3w="/home/marcgabriel/Documents/Marina_P/lncRNA_gtf_HoLDuP_run_10/contigs_sans_sens.anti.class3.quantile.0.2.gtf"

init_holdupProstate2_class3w="/home/marcgabriel/Documents/Marina_P/lncRNA_gtf_HoLDuP_run_10/contigs_sans_sens.inter.class3.quantile.0.2.gtf"


threads=8

samtools="samtools"

#path to featureCounts output
featureCounts_outputs="/home/marcgabriel/Documents/Marina_P/data_intersections/"

featureCounts_outputs="${featureCounts_outputs}/"
featureCounts_outputs=$(echo $featureCounts_outputs |sed 's/\/\//\//g')

if [ ! -d $featureCounts_outputs ]; then mkdir -p $featureCounts_outputs; fi



#### process home-made annotation ####

#convert in SAF each annotation :

#class1w prostate
final_holdupProstate1_class1w="${featureCounts_outputs}class1w_prostate1.saf"

final_holdupProstate2_class1w="${featureCounts_outputs}class1w_prostate2.saf"

#class3w prostate
final_holdupProstate1_class3w="${featureCounts_outputs}class3w_prostate1.saf"

final_holdupProstate2_class3w="${featureCounts_outputs}class3w_prostate2.saf"
 
 #class1w prostate
 paste <(awk -F'\t' '{print $9}' $init_holdupProstate1_class1w|awk -F';' '{print $1}'|awk -F' ' 'OFS="\t"{print $2}'|sed 's/\"//g'|awk 'OFS="\t"{print $1"_prostate_anti"}') <(cut -f1,4,5,7 $init_holdupProstate1_class1w)|sort -k1,1 >$final_holdupProstate1_class1w

 paste <(awk -F'\t' '{print $9}' $init_holdupProstate2_class1w|awk -F';' '{print $1}'|awk -F' ' 'OFS="\t"{print $2}'|sed 's/\"//g'|awk 'OFS="\t"{print $1"_prostate_inter"}') <(cut -f1,4,5,7 $init_holdupProstate2_class1w)|sort -k1,1 >$final_holdupProstate2_class1w
 
 #class3w prostate
 paste <(awk -F'\t' '{print $9}' $init_holdupProstate1_class3w|awk -F';' '{print $1}'|awk -F' ' 'OFS="\t"{print $2}'|sed 's/\"//g'|awk 'OFS="\t"{print $1"_prostate_anti"}') <(cut -f1,4,5,7 $init_holdupProstate1_class3w)|sort -k1,1 >$final_holdupProstate1_class3w

 paste <(awk -F'\t' '{print $9}' $init_holdupProstate2_class3w|awk -F';' '{print $1}'|awk -F' ' 'OFS="\t"{print $2}'|sed 's/\"//g'|awk 'OFS="\t"{print $1"_prostate_inter"}') <(cut -f1,4,5,7 $init_holdupProstate2_class3w)|sort -k1,1 >$final_holdupProstate2_class3w

 #remove in class 3w (prostate), IDS which are in class 1w
 
 cut -f1 $final_holdupProstate1_class1w >${featureCounts_outputs}prostate1_class3w_spe_IDs.txt
 
 cut -f1 $final_holdupProstate2_class1w >${featureCounts_outputs}prostate2_class3w_spe_IDs.txt
 

grep -v -F -f ${featureCounts_outputs}prostate1_class3w_spe_IDs.txt $final_holdupProstate1_class3w|awk 'OFS="\t"{$1=$1"_class3w";print}' >${final_holdupProstate1_class3w}.tmp && mv ${final_holdupProstate1_class3w}.tmp $final_holdupProstate1_class3w
 
grep -v -F -f ${featureCounts_outputs}prostate2_class3w_spe_IDs.txt $final_holdupProstate2_class3w|awk 'OFS="\t"{$1=$1"_class3w";print}' >${final_holdupProstate2_class3w}.tmp && mv ${final_holdupProstate2_class3w}.tmp $final_holdupProstate2_class3w


awk 'OFS="\t"{$1=$1"_class1w";print}' ${final_holdupProstate1_class1w} >${final_holdupProstate1_class1w}.tmp && mv ${final_holdupProstate1_class1w}.tmp ${final_holdupProstate1_class1w}

awk 'OFS="\t"{$1=$1"_class1w";print}' ${final_holdupProstate2_class1w} >${final_holdupProstate2_class1w}.tmp && mv ${final_holdupProstate2_class1w}.tmp ${final_holdupProstate2_class1w}

#keep only holdup transcripts with length >=200nt
awk 'OFS="\t"{mylength=($4-$3)+1;if(mylength>=200){print}}' ${final_holdupProstate1_class1w} >${final_holdupProstate1_class1w}.tmp && mv ${final_holdupProstate1_class1w}.tmp ${final_holdupProstate1_class1w}

awk 'OFS="\t"{mylength=($4-$3)+1;if(mylength>=200){print}}' ${final_holdupProstate2_class1w} >${final_holdupProstate2_class1w}.tmp && mv ${final_holdupProstate2_class1w}.tmp ${final_holdupProstate2_class1w}
 
 
awk 'OFS="\t"{mylength=($4-$3)+1;if(mylength>=200){print}}' ${final_holdupProstate1_class3w} >${final_holdupProstate1_class3w}.tmp && mv ${final_holdupProstate1_class3w}.tmp ${final_holdupProstate1_class3w}


awk 'OFS="\t"{mylength=($4-$3)+1;if(mylength>=200){print}}' ${final_holdupProstate2_class3w} >${final_holdupProstate2_class3w}.tmp && mv ${final_holdupProstate2_class3w}.tmp ${final_holdupProstate2_class3w}

#####################################
 

#### process official anotation ####

annotation_in_SAF=${featureCounts_outputs}genes.saf

if [ ! -f $annotation_in_SAF ];then

  #create SAF file (simplified annotation format) in order to give it to featureCounts
  paste <(grep "\sgene\s" $gff_gtf_ForCounting |awk 'OFS="\t"{print $9}' |awk 'OFS="\t"{split($1,a,";");for(i=1;i<=length(a);i++){if(a[i]~/^gene_id/){x=a[i]}};print x}'|awk 'OFS="\t"{sub("gene_id=","",$1);print}') <(grep "\sgene\s" $gff_gtf_ForCounting|cut -f1,4,5,7) >$annotation_in_SAF
  
  
fi


#### merge official annotation with holdup annotation  ####

#gff format
#give a gene_type to the home-made annotation, and merge them with the official annotation in gff format
 cat <(cat $final_holdupProstate1_class1w $final_holdupProstate2_class1w $final_holdupProstate1_class3w $final_holdupProstate2_class3w|awk 'OFS="\t"{print $2,"holdup","gene",$3,$4,".",$5,".","ID="$1";gene_id="$1";gene_name="$1";"}') <(grep -v "^#" $gff_gtf_ForCounting|grep -E "gene\s")|awk 'OFS="\t"{if($9~"_anti_class3w"){$9=$9"gene_type=lncRNA_anti_class3w"}else if($9~"_inter_class3w"){$9=$9"gene_type=lncRNA_inter_class3w"}else if($9~"_anti_class1w"){$9=$9"gene_type=lncRNA_anti_class1w"}else if($9~"_inter_class1w"){$9=$9"gene_type=lncRNA_inter_class1w"};print}' >${featureCounts_outputs}extended_annotation.gff
 
 #ID=
 awk '{print $9}' ${featureCounts_outputs}extended_annotation.gff|awk 'OFS="\t"{split($1,a,";");split($1,b,";");for(i=1;i<=length(a);i++){if(a[i]~/^gene_type/){x=a[i]};if(a[i]~/^gene_id=/){y=b[i]}};print y,x}'|awk 'OFS="\t"{sub("gene_id=","",$1);sub("gene_type=","",$2);print}' >${featureCounts_outputs}matching_genes_types.tsv
 
 #link each gene ID to a defined type
 join -t $'\t' -12 -21 <(sort -k2,2 ${featureCounts_outputs}matching_genes_types.tsv) <(sort -k1,1 $gene_type)|cut -f 2,3|awk 'OFS="\t"{if($1~"_prostate_"){$3="prostate"}else if($1~"ENSG00"){$3="gencode"};print}' >${featureCounts_outputs}matching_genes_types.tmp && mv ${featureCounts_outputs}matching_genes_types.tmp ${featureCounts_outputs}matching_genes_types.tsv
 
 
 sed -i '1 i\'"ID""\t""type""\t""experiment"'' ${featureCounts_outputs}matching_genes_types.tsv
 
 #in saf format
 cat $final_holdupProstate1_class1w $final_holdupProstate2_class1w $final_holdupProstate1_class3w $final_holdupProstate2_class3w $annotation_in_SAF |sort -k1,1 |cat <(echo -e "GeneID\tChr\tStart\tEnd\tStrand") - >${featureCounts_outputs}extended_annotation.saf
 
 
 grep -E "ENSG" ${featureCounts_outputs}extended_annotation.gff >${featureCounts_outputs}gencode_only.gff
 
 grep -E "\sholdup\s" ${featureCounts_outputs}extended_annotation.gff|grep "class1w" >${featureCounts_outputs}holdup_prostate_class1w_only.gff
 
 grep -E "\sholdup\s" ${featureCounts_outputs}extended_annotation.gff|grep "class3w" >${featureCounts_outputs}holdup_prostate_class3w_only.gff
 
 grep -E "\sholdup\s" ${featureCounts_outputs}extended_annotation.gff|grep "class1w" |grep "_inter_">${featureCounts_outputs}holdup_prostate_class1w_only_inter_only.gff
 
 grep -E "\sholdup\s" ${featureCounts_outputs}extended_annotation.gff|grep "class1w" |grep "_anti_">${featureCounts_outputs}holdup_prostate_class1w_only_anti_only.gff
 
 grep -E "\sholdup\s" ${featureCounts_outputs}extended_annotation.gff|grep "class3w"|grep "_inter_" >${featureCounts_outputs}holdup_prostate_class3w_only_inter_only.gff
 
 grep -E "\sholdup\s" ${featureCounts_outputs}extended_annotation.gff|grep "class3w"|grep "_anti_" >${featureCounts_outputs}holdup_prostate_class3w_only_anti_only.gff
 
#exit

##########################################

##########################################

#normal files
normal_files=($(find "/media/marcgabriel/SAMSUNG/24_prostate_files_CURIE/" -name "*unique.bam"|grep -E "B65T2\.|B65T4\.|B65T6\.|B65T8\.|B67T10\.|B67T13\.|B67T14\.|B67T15\.")) || { echo "no files (normal) !!" 1>&2; exit; }

#tumoral files
tumoral_files=($(find "/media/marcgabriel/SAMSUNG/24_prostate_files_CURIE/" -name "*unique.bam"|grep -E "B67T18\.|B67T19\.|B67T20\.|B67T22\.|B65T1\.|B65T3\.|B65T5\.|B65T7\.|B67T11\.|B67T17\.|B67T21\.|B67T12\.|B67T23\.|B67T24\.|B67T9\.|B67T16\.")) || { echo "no files (tumoral) !!" 1>&2; exit; }

echo -e "\nnb normal files : ${#normal_files[*]}"
  
echo -e "\nnb tumoral files : ${#tumoral_files[*]}"

all_files=(${normal_files[*]} ${tumoral_files[*]})

echo -e "\ntotal files to process : ${#all_files[*]}"


#remove existing summary file
if [ -f ${featureCounts_outputs}counts_summary.txt ];then rm -rf ${featureCounts_outputs}counts_summary.txt;fi

for i in ${all_files[*]};do

  files_name=$(basename $i|sed 's/\.mapped\.unique\.bam//g')
  
  
  echo -e "\n-- sample : ${files_name}, file $i --\n"

  #genes
  if [ ! -f ${featureCounts_outputs}${files_name}_counts.tsv ];then
  
	  #determine reads number, in order to compute the RPKM/FPKM
	  mapped_paired=$($samtools view -F 0x8 -F 0x4 -f 0x1 -c $i)
	  
	  mapped_singletons=$($samtools view -f 0x8 -c $i)
	  
	  total_reads=$(awk -v mapped_paired=$mapped_paired -v mapped_singletons=$mapped_singletons 'BEGIN{paired=mapped_paired/2;print paired+mapped_singletons}')
	  
	  echo -e "\n\t- total fragments is : $total_reads"
	  
	  echo -e "${files_name}\t${total_reads}" >>${featureCounts_outputs}counts_summary.txt
	  
	  start=$(date)
	  
	  #-O allowMultiOverlap (reads on multiple features are counted for each)
	  $featureCounts_prog -F "SAF" -p $orientation -T $threads -O -a ${featureCounts_outputs}extended_annotation.saf -o ${featureCounts_outputs}${files_name}_counts.tsv $i
	  
	  end=$(date)
	  
	  echo -e "\n\t- time for ${files_name} : $start - $end\n"
	  
	  #$7+1
	  grep -v "^#" ${featureCounts_outputs}${files_name}_counts.tsv|awk -v total_reads=$total_reads 'NR>1{OFS="\t";RPKM=((($7)*1000*1000000)/($6*total_reads));print $0,RPKM}' |tee ${featureCounts_outputs}${files_name}_RPKM.tsv|cut -f1,7 >${featureCounts_outputs}${files_name}_counts.tmp && mv ${featureCounts_outputs}${files_name}_counts.tmp ${featureCounts_outputs}${files_name}_counts.tsv
	  
	  
  fi

  
done







