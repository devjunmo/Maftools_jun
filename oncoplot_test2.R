library(maftools)
library(tibble)

#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

tibble(read.delim(laml.maf))


clin_test <- tibble(read.delim(laml.clin))
clin_test





laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
laml

#Shows sample summry.
getSampleSummary(laml)
#Shows gene summary.
getGeneSummary(laml)
#shows clinical data associated with samples
getClinicalData(laml)
#Shows all fields in MAF
getFields(laml)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')

plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)


# 7.2 Oncoplots


#oncoplot for top ten mutated genes.
oncoplot(maf = laml, top = 10)

# 7.3 Transition and Transversions.
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)


# 7.4 Lollipop plots for amino acid changes

#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
lollipopPlot(
  maf = laml,
  gene = 'DNMT3A',
  AACol = 'Protein_Change',
  showMutationRate = TRUE,
  labelPos = 882
)

?plotProtein
# plotProtein(gene = "TP53", refSeqID = "NM_000546")


# 7.5 Rainfall plots

brca <- system.file("extdata", "brca.maf.gz", package = "maftools")
brca = read.maf(maf = brca, verbose = FALSE)


rainfallPlot(maf = brca, detectChangePoints = TRUE, pointSize = 0.4)


# 7.6 Compare mutation load against TCGA cohorts

aml.mutload = tcgaCompare(maf = laml, cohortName = 'Example-LAML', logscale = TRUE, capture_size = 50)

# 7.7 Plotting VAF

plotVaf(maf = laml, vafCol = 'i_TumorVAF_WU')

# 8 Processing copy-number data

all.lesions <- system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes <- system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
scores.gis <- system.file("extdata", "scores.gistic", package = "maftools")

laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)

laml.gistic

gisticChromPlot(gistic = laml.gistic, markBands = "all")


# 8.2.2 Bubble plot

gisticBubblePlot(gistic = laml.gistic)

# 8.2.3 oncoplot

gisticOncoPlot(gistic = laml.gistic, clinicalData = getClinicalData(x = laml), clinicalFeatures = 'FAB_classification', sortByAnnotation = TRUE, top = 10)

# 8.2.4 Visualizing CBS segments

tcga.ab.009.seg <- system.file("extdata", "TCGA.AB.3009.hg19.seg.txt", package = "maftools")
plotCBSsegments(cbsFile = tcga.ab.009.seg)

# 9 Analysis

#exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))

# 9.2 Detecting cancer driver genes based on positional clustering

laml.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')

head(laml.sig)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

# 9.3 Adding and summarizing pfam domains

laml.pfam = pfamDomains(maf = laml, AACol = 'Protein_Change', top = 10)
laml.pfam$proteinSummary[,1:7, with = FALSE]
laml.pfam$domainSummary[,1:3, with = FALSE]
#

# 9.4 Survival analysis
# 9.4.1 Mutation in any given genes

#Survival analysis based on grouping of DNMT3A mutation status
getClinicalData(laml)
mafSurvival(maf = laml, genes = 'DNMT3A', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)


# 9.4.2 Predict genesets associated with survival
#Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
prog_geneset = survGroup(maf = laml, top = 20, geneSetSize = 2, time = "days_to_last_followup", Status = "Overall_Survival_Status", verbose = FALSE)

print(prog_geneset)
mafSurvGroup(maf = laml, geneSet = c("DNMT3A", "FLT3"), time = "days_to_last_followup", Status = "Overall_Survival_Status")



# 9.5 Comparing two cohorts (MAFs) # 

#Primary APL MAF
primary.apl = system.file("extdata", "APL_primary.maf.gz", package = "maftools")
primary.apl
tibble(read.delim(primary.apl))

primary.apl = read.maf(maf = primary.apl)
primary.apl = read.maf(maf = tibble(read.delim(primary.apl)))

#Relapse APL MAF
relapse.apl = system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
tibble(read.delim(relapse.apl))
relapse.apl = read.maf(maf = relapse.apl)
relapse.apl = read.maf(maf = tibble(read.delim(relapse.apl)))

pt.vs.rt <- mafCompare(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary', m2Name = 'Relapse', minMut = 5)
print(pt.vs.rt)

forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1)


# 9.5.2 Co-onco plots

?coOncoplot
genes = c("PML", "RARA", "RUNX1", "ARID1B", "FLT3") # 위에서 코호트간 뮤테이션 카운트 차이가 컸던 상위 5개 genes
coOncoplot(m1 = primary.apl, m2 = relapse.apl, m1Name = 'PrimaryAPL', m2Name = 'RelapseAPL', genes = genes, removeNonMutated = TRUE)


# 9.5.3 Co-bar plots

coBarplot(m1 = primary.apl, m2 = relapse.apl, m1Name = "Primary", m2Name = "Relapse")



# 9.5.4 Lollipop plot-2

lollipopPlot2(m1 = primary.apl, m2 = relapse.apl, gene = "PML", AACol1 = "amino_acid_change", AACol2 = "amino_acid_change", m1_name = "Primary", m2_name = "Relapse")

# 9.6 Clinical enrichment analysis

fab.ce = clinicalEnrichment(maf = laml, clinicalFeature = 'FAB_classification')
fab.ce$groupwise_comparision[p_value < 0.05]
plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)

# 9.7 Drug-Gene Interactions
dgi = drugInteractions(maf = laml, fontSize = 0.75)

dnmt3a.dgi = drugInteractions(genes = "DNMT3A", drugs = TRUE)

#Printing selected columns.
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]




# 9.8 Oncogenic Signaling Pathways

OncogenicPathways(maf = laml)

PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS")

























































