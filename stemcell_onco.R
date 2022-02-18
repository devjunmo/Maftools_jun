library(maftools)
library(tibble)
library(dplyr)

# root_dir = 'E:/stemcell/somatic_analysis/maf/mutect2/DP_AF_filtered_maf/exonic_maf/cohort_comp'
# root_dir = 'E:/stemcell/somatic_analysis/maf/mutect2/DP_AF_filtered_maf/exonic_maf/diff_comp/cohort_comp'
root_dir = 'E:/stemcell/somatic_analysis/maf/mutect2/DP_AF_filtered_maf/Diff_whole_regions/cohort_comp'
# maf_f_name = 'ips_cohort_maf_for_cooncoplot.maf'
maf_f_name = 'diff_cohort_for_oncoplot.maf'
# clin_f_name = 'IPS_cohort_info.csv'
clin_f_name = 'diff_cohort_info.csv'

getwd()
old.path <- setwd(root_dir)
getwd()

maf_df <-  tibble(read.delim(maf_f_name, na.strings = c("")))
clin_df <- tibble(read.csv(clin_f_name, fileEncoding = 'UTF-8-BOM',
                           na.strings = c("")))

print(maf_df, n=3, width = Inf)
print(clin_df, n=3, width = Inf)
clin_df

laml = read.maf(maf = maf_df,
                clinicalData = clin_df,
                verbose = FALSE)

?oncoplot


draw_oncoplot <- function(maf_obj_, 
                          clin_feature_vec_,
                          annotationOrder_vec_,
                          title_txt_,
                          top_genes_count_)
  {
  
  oncoplot(maf = maf_obj_, top =top_genes_count_, draw_titv = F, drawRowBar=F,
           clinicalFeatures=clin_feature_vec_,
           showTumorSampleBarcodes=T, barcode_mar=7, barcodeSrt=45, gene_mar=10,
           annotationOrder=annotationOrder_vec_,
           legend_height=4, sortByAnnotation=T, removeNonMutated=F, drawBox=T,
           titleText=title_txt_, titleFontSize=1.6,
           anno_height=0.8, SampleNamefontSize=1.2)
}



# Origin별 소팅
# oncoplot(maf = laml, top =50, draw_titv = F, drawRowBar=F,
#          clinicalFeatures=c('Origin_sample_Barcode', 'Donor_age', 'Mt_mutation'),
#          showTumorSampleBarcodes=T, barcode_mar=5, barcodeSrt=30, gene_mar=10,
#          annotationOrder=c('AM168', 'AM172', 'AM179', 'Blo86', 'Blo172',
#                            'Fib86', 'Fib172'),
#          legend_height=4, sortByAnnotation=T, removeNonMutated=F, drawBox=T,
#          titleText='IPS Oncoplot: sort by origin samples', titleFontSize=1.6,
#          anno_height=0.8, SampleNamefontSize=1.2)


# IPS
clinicalFeatures=c('Origin_sample_Barcode', 'Donor_age', 'Mt_mutation')
annotationOrder=c('AM168', 'AM172', 'AM179', 'Blo86', 'Blo172',
                  'Fib86', 'Fib172')

titleText='IPS Oncoplot: sort by origin samples'

draw_oncoplot(laml, clinicalFeatures, annotationOrder, titleText)


# Tera
clinicalFeatures=c('Origin_sample', 'Donor_age', 'Mt_mutation')
annotationOrder=c('AM168', 'AM172', 'AM179', 'Blo86', 'Blo172',
                  'Fib86')
titleText='Tera Oncoplot: sort by origin samples'

draw_oncoplot(laml, clinicalFeatures, annotationOrder, titleText, 50)


# Diff
clinicalFeatures=c('Origin', 'Mt_mutation')
annotationOrder=c('Blo172', 'hiPS29-A-p49-1', 'hiPS29-B-p49-1')
# titleText='Oncoplot: Differentiated cell samples'
titleText='Oncoplot: Differentiated cell samples (Whole regions)'

draw_oncoplot(laml, clinicalFeatures, annotationOrder, titleText, 50)


# Diff-rm_tera

maf_df <-  tibble(read.delim(maf_f_name, na.strings = c("")))
unique(maf_df$Tumor_Sample_Barcode)

library(dplyr)

maf_df_rmTera <- filter(maf_df, Tumor_Sample_Barcode!="29-A-P50-Tera" &
                          Tumor_Sample_Barcode!="29-B-P50-Tera")
maf_df_rmTera[maf_df_rmTera$Tumor_Sample_Barcode=="29-A-Beta-cell",]

maf_df_rmTera
maf_df_rmTera[maf_df_rmTera$Tumor_Sample_Barcode=="29-B-Beta-cell", ]

laml.rmTera = read.maf(maf = maf_df_rmTera,
                clinicalData = clin_df,
                verbose = FALSE)

clinicalFeatures=c('Origin', 'Mt_mutation')
annotationOrder=c('Blo172', 'hiPS29-A-p49-1', 'hiPS29-B-p49-1')
titleText='Oncoplot: Differentiated cell samples (Whole regions)'

draw_oncoplot(laml.rmTera, clinicalFeatures, annotationOrder, titleText, 50)



# only snp

# maf_df_rmTera_SNP <- filter(maf_df, Tumor_Sample_Barcode!="29-A-P50-Tera" &
#                           Tumor_Sample_Barcode!="29-B-P50-Tera" &
#                             Tumor_Sample_Barcode!="hiPS29-A-p49-1" &
#                             Tumor_Sample_Barcode!="hiPS29-B-p49-1" &
#                             Tumor_Sample_Barcode!="hiPS29A-e-p-DNA" &
#                             Tumor_Sample_Barcode!="hiPS29B-e-p-DNA" &
#                             Variant_Type=='SNP')

maf_df_rmTera_SNP <- filter(maf_df, Tumor_Sample_Barcode!="29-A-P50-Tera" &
                              Tumor_Sample_Barcode!="29-B-P50-Tera" &
                              Variant_Type=='SNP')

maf_df_rmTera_SNP

laml.rmTera.SNP = read.maf(maf = maf_df_rmTera_SNP,
                       clinicalData = clin_df,
                       verbose = FALSE)

clinicalFeatures=c('Origin', 'Mt_mutation')
annotationOrder=c('Blo172', 'hiPS29-A-p49-1', 'hiPS29-B-p49-1')
titleText='Oncoplot: Differentiated cell samples & only SNP (Whole regions)'

draw_oncoplot(laml.rmTera.SNP, clinicalFeatures, annotationOrder, titleText, 50)




setwd(old.path)
getwd()























