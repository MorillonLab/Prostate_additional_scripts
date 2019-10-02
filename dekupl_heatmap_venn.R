#!/usr/bin/env Rscript

library(ComplexHeatmap)

library(gplots)

library(circlize)


suppressMessages(library(plotrix))

suppressMessages(library(Biostrings))


### input data ###

home<-"/home/marcgabriel/Documents/Marina_P/"

complete_table<-"/home/marcgabriel/Documents/Marina_P/DiffContigsInfos.tsv"

complete_table<-read.delim(file=complete_table,sep = "\t",check.names = F,header=T)

dekupl_contigs<-list("/home/marcgabriel/Documents/Marina_P/kmers_annotations/kmers_contiguous.gff",
"/home/marcgabriel/Documents/Marina_P/kmers_annotations/kmers_spliced.gff",
"/home/marcgabriel/Documents/Marina_P/kmers_annotations/kmers_repeat.gff",
"/home/marcgabriel/Documents/Marina_P/kmers_annotations/kmers_unmapped.gff")

#################

#"kmer_81545"
probes<-c("ctg_28650",
          "ctg_61528",
          "ctg_61472",
          "ctg_117356",
          "ctg_9446",
          "ctg_44030",
          "ctg_512",
          "ctg_57223",
          "ctg_17297",
          "ctg_111348",
          "ctg_63866" ,
          "ctg_29077",
          "ctg_2815",
          "ctg_119680",
          "ctg_36195",
          "ctg_123090",
          "ctg_73782",
          "ctg_25348",
          "ctg_111158",
          "ctg_105149",
          "ctg_23999",
          "ctg_37852",
          "ctg_104447")

subset_table<-data.frame()
for(i in 1:length(dekupl_contigs)){
  
  mytable<-read.delim(file=dekupl_contigs[[i]],sep = "\t",check.names = F,header=F)
  
  mytable<-mytable[c(1,4,5,7,9)]
  
  names(mytable)<-c("chromosome","start","end","strand","ID")
  
  subset_table<-rbind(subset_table,mytable)
  
}

subset_table<-merge(subset_table["ID"],complete_table[,!names(complete_table)%in%"ID"], 
                    by.x="ID",
                    by.y="LineInSam",
                    all.x=TRUE,
                    all.y=FALSE)



samples_col<-24

cond1_col<-8

start_samples<-(ncol(subset_table)-samples_col)+1

end_cond1<-(start_samples+cond1_col)-1

start_cond2<-end_cond1+1

names(subset_table)[start_samples:end_cond1]<-paste(rep("normal",cond1_col),seq(1:cond1_col),sep="_")

names(subset_table)[start_cond2:ncol(subset_table)]<-paste(rep("tumoral",samples_col-cond1_col),seq(1:(samples_col-cond1_col)),sep="_")

my_palette<-c("darkblue","blue","cornflowerblue","lightblue","white","salmon","#E44D2E","red","darkred")


my_matrix<-as.matrix(subset_table[,-c(1:start_samples-1)])

min_val<-min(my_matrix[my_matrix>0])/2

my_matrix[my_matrix==0]<-min_val

my_matrix<-log10(my_matrix)

#ha <-HeatmapAnnotation(cn = anno_text(colnames(my_matrix), just = "top", offset = unit(1, "npc"),rot=20),which ="column")

#ha_height = max_text_height(colnames(my_matrix)) 

hb<-columnAnnotation(df =data.frame(conditions = c(rep("normal",cond1_col),rep("tumoral",samples_col-cond1_col))),
                      col =list(conditions =c("normal" = "grey60", "tumoral" = "purple")),annotation_legend_param = list(labels_gp=gpar(fontsize=16),title_gp=gpar(fontsize = 16),grid_width= unit(10, "mm"),nrow = 1, by_row = TRUE,title_position = "topcenter"))

rownames(my_matrix)<-paste("ctg",subset_table$ID,sep="_")

subset = which(rownames(my_matrix)%in%probes, arr.ind=T)

labels = rownames(my_matrix)[subset]

hc<-rowAnnotation(link = row_anno_link(at = subset, labels = labels),
              width = unit(1, "cm") + max_text_width(labels))

#png(filename=paste(home,"heatmap_dekupl_contigs_expression.png",sep=""),width=1000,height=800)

  ht1 = Heatmap(my_matrix, col=colorRamp2(c(seq(min(my_matrix), 1, len=5),1.5,2,2.5,3),my_palette),
                name = "log10(normalized expression)",
                column_dend_side = "bottom",
                column_title = "", bottom_annotation = hb,
                row_dend_width = unit(3, "cm"),
                clustering_method_rows="ward.D",
                clustering_method_columns="ward.D",
                clustering_distance_rows = "euclidean",clustering_distance_columns = "euclidean",
                cluster_columns =T,cluster_rows =T,show_row_names=F,show_column_names = F,
                heatmap_legend_param=list(labels_gp=gpar(fontsize=12),title_gp=gpar(fontsize = 14),grid_width= unit(5, "mm")))
  hm_list<-ht1+hc
  draw(hm_list, heatmap_legend_side = "right",row_title=paste(nrow(my_matrix)," Dekupl contigs",sep=""),annotation_legend_side = "bottom")

#dev.off()

#matrix for just the selected contigs/probes
my_matrix<-my_matrix[which(row.names(my_matrix)%in%probes),]

ht1 = Heatmap(my_matrix, col=colorRamp2(c(seq(min(my_matrix), 1, len=5),1.5,2,2.5,3),my_palette),
              name = "log10(normalized expression)",
              column_dend_side = "bottom",
              column_title = "", bottom_annotation = hb,
              row_dend_width = unit(3, "cm"),
              clustering_method_rows="ward.D",
              clustering_method_columns="ward.D",
              clustering_distance_rows = "euclidean",clustering_distance_columns = "euclidean",
              row_names_gp = gpar(fontsize = 20),
              cluster_columns =T,cluster_rows =T,show_row_names=T,show_column_names = F,
              heatmap_legend_param=list(labels_gp=gpar(fontsize=16),
                                        title_gp=gpar(fontsize = 16),grid_width= unit(10, "mm"))
              )



hm_list<-ht1

png(filename=paste(home,"heatmap_",nrow(my_matrix),"_dekupl_contigs_expression.png",sep=""),width=800,height=1000)

draw(hm_list, heatmap_legend_side = "right",row_title=paste(nrow(my_matrix)," Dekupl contigs",sep=""),annotation_legend_side = "bottom")
dev.off()

pdf(paste(home,"heatmap_",nrow(my_matrix),"_dekupl_contigs_expression.pdf",sep=""),width=10,height=10)
draw(hm_list, heatmap_legend_side = "right",row_title=paste(nrow(my_matrix)," Dekupl contigs",sep=""),annotation_legend_side = "bottom")
dev.off()

#postscript(paste(home,"heatmap_",nrow(my_matrix),"_dekupl_contigs_expression.eps",sep=""), horizontal = FALSE, onefile = FALSE, paper = "special")
cairo_ps(paste(home,"heatmap_",nrow(my_matrix),"_dekupl_contigs_expression.eps",sep=""),width=10,height=10)

 draw(hm_list, heatmap_legend_side = "right",row_title=paste(nrow(my_matrix)," Dekupl contigs",sep=""),annotation_legend_side = "bottom")

dev.off()


############

nb_exonic_contigs<-nrow(subset_table[which(subset_table$exonic==T & subset_table$UTR==F),])

mapped_IDs<-subset_table[which(subset_table$exonic==T & subset_table$UTR==F),]$ID

nb_intronic_contigs<-nrow(subset_table[which(subset_table$exonic==F & subset_table$intronic==T & subset_table$UTR==F),])

mapped_IDs<-c(mapped_IDs,subset_table[which(subset_table$exonic==F & subset_table$intronic==T & subset_table$UTR==F),]$ID)

nb_UTR_contigs<-nrow(subset_table[which(subset_table$UTR==T),])

mapped_IDs<-c(mapped_IDs,subset_table[which(subset_table$UTR==T),]$ID)

antisense_contigs<-subset_table[which(subset_table$HUGO_ID_as_gene!="none"),]

nb_antisense_contigs<-nrow(antisense_contigs[which(!antisense_contigs$ID%in%mapped_IDs),])

mapped_IDs<-c(mapped_IDs,antisense_contigs[which(!antisense_contigs$ID%in%mapped_IDs),]$ID)

nb_intergenic_contigs<-nrow(subset_table[which(!subset_table$ID%in%mapped_IDs & subset_table$is_mapped==T),])


contigs_classification<-data.frame(features=c(paste("exonic\n",round((nb_exonic_contigs/nrow(subset_table))*100,2),"%",sep=""),
                                              paste("intronic\n",round((nb_intronic_contigs/nrow(subset_table))*100,2),"%",sep=""),
                                              paste("UTR\n",round((nb_UTR_contigs/nrow(subset_table))*100,2),"%",sep=""),
                                              paste("intergenic\n",round((nb_intergenic_contigs/nrow(subset_table))*100,2),"%",sep=""),
                                              paste("antisense\n",round((nb_antisense_contigs/nrow(subset_table))*100,2),"%",sep="")),
                                   values=c(round((nb_exonic_contigs/nrow(subset_table))*100,2),
                                            round((nb_intronic_contigs/nrow(subset_table))*100,2),
                                            round((nb_UTR_contigs/nrow(subset_table))*100,2),
                                            round((nb_intergenic_contigs/nrow(subset_table))*100,2),
                                            round((nb_antisense_contigs/nrow(subset_table))*100,2)),
                                   color=c("cornflowerblue","red","purple","grey60","darkgreen"))


png(filename=paste(home,"piechart_dekupl_contigs_distribution.png",sep=""),width=1000,height=800)

  pie(contigs_classification$values,contigs_classification$features,col=as.character(contigs_classification$color),
      main=paste("Dekupl contigs distribution across genomic features\nTotal contigs = ",format(nrow(subset_table),big.mark=" "),"\nAllocation priority : UTR > exonic > intronic > antisense > intergenic",sep=""))

dev.off()

