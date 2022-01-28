
# oncoplot

library(maftools)

#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')

test_clin <- read.delim(laml.clin)
test_clin

laml = read.maf(maf = laml.maf,
                clinicalData = laml.clin,
                verbose = FALSE)

# 0.1 Including Transition/Transversions into oncoplot

oncoplot(maf = laml, draw_titv = TRUE)


# 0.2 Changing colors for variant classifications

#One can use any colors, here in this example color palette from RColorBrewer package is used
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

print(vc_cols)
#>   Frame_Shift_Del Missense_Mutation Nonsense_Mutation         Multi_Hit 
#>         "#A6CEE3"         "#1F78B4"         "#B2DF8A"         "#33A02C" 
#>   Frame_Shift_Ins      In_Frame_Ins       Splice_Site      In_Frame_Del 
#>         "#FB9A99"         "#E31A1C"         "#FDBF6F"         "#FF7F00"

oncoplot(maf = laml, colors = vc_cols, top = 10)

# 0.3 Including copy number data into oncoplots.

#GISTIC results LAML
all.lesions =
  system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes =
  system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes =
  system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
scores.gis =
  system.file("extdata", "scores.gistic", package = "maftools")

#Read GISTIC results along with MAF
laml.plus.gistic = read.maf(
  maf = laml.maf,
  gisticAllLesionsFile = all.lesions,
  gisticAmpGenesFile = amp.genes,
  gisticDelGenesFile = del.genes,
  gisticScoresFile = scores.gis,
  isTCGA = TRUE,
  verbose = FALSE, 
  clinicalData = laml.clin
)

# This plot shows frequent deletions in TP53 gene which is located on one of the significantly deleted locus 17p13.2.

oncoplot(maf = laml.plus.gistic, top = 10)



# 0.3.2 Custom copy-number table

set.seed(seed = 1024)
barcodes = as.character(getSampleSummary(x = laml)[,Tumor_Sample_Barcode])
#Random 20 samples
dummy.samples = sample(x = barcodes,
                       size = 20,
                       replace = FALSE)

#Genarate random CN status for above samples
cn.status = sample(
  x = c('ShallowAmp', 'DeepDel', 'Del', 'Amp'),
  size = length(dummy.samples),
  replace = TRUE
)

custom.cn.data = data.frame(
  Gene = "DNMT3A",
  Sample_name = dummy.samples,
  CN = cn.status,
  stringsAsFactors = FALSE
)

head(custom.cn.data)
#>     Gene  Sample_name         CN
#> 1 DNMT3A TCGA-AB-2898 ShallowAmp
#> 2 DNMT3A TCGA-AB-2879        Del
#> 3 DNMT3A TCGA-AB-2920        Amp
#> 4 DNMT3A TCGA-AB-2866        Del
#> 5 DNMT3A TCGA-AB-2892        Del
#> 6 DNMT3A TCGA-AB-2863 ShallowAmp

laml.plus.cn = read.maf(maf = laml.maf,
                        cnTable = custom.cn.data,
                        verbose = FALSE)

oncoplot(maf = laml.plus.cn, top = 5)



# 0.4 Bar plots

#Selected AML driver genes
aml_genes = c("TP53", "WT1", "PHF6", "DNMT3A", "DNMT3B", "TET1", "TET2", "IDH1", "IDH2", "FLT3", "KIT", "KRAS", "NRAS", "RUNX1", "CEBPA", "ASXL1", "EZH2", "KDM6A")

#Variant allele frequcnies (Right bar plot)
aml_genes_vaf = subsetMaf(maf = laml, genes = aml_genes, fields = "i_TumorVAF_WU", mafObj = FALSE)[,mean(i_TumorVAF_WU, na.rm = TRUE), Hugo_Symbol]
colnames(aml_genes_vaf)[2] = "VAF"
head(aml_genes_vaf)
#>    Hugo_Symbol      VAF
#> 1:       ASXL1 37.11250
#> 2:       CEBPA 22.00235
#> 3:      DNMT3A 43.51556
#> 4:      DNMT3B 37.14000
#> 5:        EZH2 68.88500
#> 6:        FLT3 34.60294

#MutSig results (Right bar plot)
laml.mutsig = system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
laml.mutsig = data.table::fread(input = laml.mutsig)[,.(gene, q)]
laml.mutsig[,q := -log10(q)] #transoform to log10
head(laml.mutsig)
#>      gene        q
#> 1:   FLT3 12.64176
#> 2: DNMT3A 12.64176
#> 3:   NPM1 12.64176
#> 4:   IDH2 12.64176
#> 5:   IDH1 12.64176
#> 6:   TET2 12.64176


oncoplot(
  maf = laml,
  genes = aml_genes,
  leftBarData = aml_genes_vaf,
  leftBarLims = c(0, 100),
  rightBarData = laml.mutsig,
  rightBarLims = c(0, 20)
)


# 0.5 Including annotations

getClinicalData(x = laml)
#>      Tumor_Sample_Barcode FAB_classification days_to_last_followup
#>   1:         TCGA-AB-2802                 M4                   365
#>   2:         TCGA-AB-2803                 M3                   792
#>   3:         TCGA-AB-2804                 M3                  2557
#>   4:         TCGA-AB-2805                 M0                   577
#>   5:         TCGA-AB-2806                 M1                   945
#>  ---                                                              
#> 189:         TCGA-AB-3007                 M3                  1581
#> 190:         TCGA-AB-3008                 M1                   822
#> 191:         TCGA-AB-3009                 M4                   577
#> 192:         TCGA-AB-3011                 M1                  1885
#> 193:         TCGA-AB-3012                 M3                  1887
#>      Overall_Survival_Status
#>   1:                       1
#>   2:                       1
#>   3:                       0
#>   4:                       1
#>   5:                       1
#>  ---                        
#> 189:                       0
#> 190:                       1
#> 191:                       1
#> 192:                       0
#> 193:                       0

oncoplot(maf = laml, genes = aml_genes, clinicalFeatures = 'FAB_classification')

#Color coding for FAB classification
fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
names(fabcolors) = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7")
fabcolors = list(FAB_classification = fabcolors)

print(fabcolors)
#> $FAB_classification
#>        M0        M1        M2        M3        M4        M5        M6        M7 
#> "#D53E4F" "#F46D43" "#FDAE61" "#FEE08B" "#E6F598" "#ABDDA4" "#66C2A5" "#3288BD"

oncoplot(
  maf = laml, genes = aml_genes,
  clinicalFeatures = 'FAB_classification',
  sortByAnnotation = TRUE,
  annotationColor = fabcolors
)

getFields(x = laml)

oncoplot(maf = laml, genes = aml_genes,
         additionalFeature = c("Tumor_Seq_Allele2", "C"))



# 0.7 Group by Pathways

oncoplot(maf = laml, pathways = "auto", gene_mar = 8, fontSize = 0.6)



# 0.7.2 Custom pathways

pathways = data.frame(
  Genes = c(
    "TP53",
    "WT1",
    "PHF6",
    "DNMT3A",
    "DNMT3B",
    "TET1",
    "TET2",
    "IDH1",
    "IDH2",
    "FLT3",
    "KIT",
    "KRAS",
    "NRAS",
    "RUNX1",
    "CEBPA",
    "ASXL1",
    "EZH2",
    "KDM6A"
  ),
  Pathway = rep(c(
    "TSG", "DNAm", "Signalling", "TFs", "ChromMod"
  ), c(3, 6, 4, 2, 3)),
  stringsAsFactors = FALSE
)
pathways
head(pathways)

oncoplot(maf = laml, pathways = pathways, gene_mar = 8, fontSize = 0.6)



# 0.8 Combining everything

oncoplot(
  maf = laml.plus.gistic,
  draw_titv = TRUE,
  pathways = pathways,
  clinicalFeatures = c('FAB_classification', 'Overall_Survival_Status'),
  sortByAnnotation = TRUE,
  additionalFeature = c("Tumor_Seq_Allele2", "C"),
  leftBarData = aml_genes_vaf,
  leftBarLims = c(0, 100),
  rightBarData = laml.mutsig[,.(gene, q)],
)
























