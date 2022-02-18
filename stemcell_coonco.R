  #############################################
 # stemcell: processing mafs for co-oncoplot #
#############################################


# maf에 코호트 정보 추가
# 필요한 컬럼만 추출

## -- python -- (VCF_with_py에서 진행)

library(maftools)
library(tibble)
library(dplyr)

# root_dir = 'E:/stemcell/somatic_analysis/maf/mutect2/DP_AF_filtered_maf/exonic_maf/cohort_comp'
root_dir = 'E:/stemcell/somatic_analysis/maf/mutect2/DP_AF_filtered_maf/exonic_maf/diff_comp/cohort_comp'
maf_f_name = 'diff_cohort_for_oncoplot.maf'
old.path <- setwd(root_dir)
getwd()

maf_df <-  tibble(read.delim(maf_f_name))
# maf_df <- read.delim(maf_f_name)

maf_df

?coOncoplot

draw_co_oncoplot <- function(df1, df2, name1, name2){
  
  # comp_obj <- mafCompare(m1 = primary.apl, m2 = relapse.apl, 
  #                        m1Name = 'Primary', m2Name = 'Relapse', minMut = 5)
  # print(comp_obj)
  
  maf1 = read.maf(maf = df1)
  maf2 = read.maf(maf = df2)
  
  
  coOncoplot(m1 = maf1, m2 = maf2, m1Name = name1, 
             m2Name = name2, removeNonMutated = F,
             showSampleNames=TRUE, SampleNamefont=0.9,
             legend_height=1.3, barcode_mar=8, outer_mar=3)
}



# Diff sample 단위 비교

## IPS29A 유래 vs IPS29B 유래

draw_co_oncoplot(maf_df[maf_df$Origin== 'hiPS29-A-p49-1', ],
                 maf_df[maf_df$Origin == 'hiPS29-B-p49-1', ],
                 'From_IPS29A', 'From_IPS29B')
# ?mafCompare




# 코호트간 비교 
# 컬럼명, 대상값1, 대상값2 -> 
# ?coOncoplot

# origin sample 단위 비교

## AM168 vs AM172

draw_co_oncoplot(maf_df[maf_df$Origin_sample_Barcode == 'AM168', ],
                 maf_df[maf_df$Origin_sample_Barcode == 'AM172', ],
                 'AM168', 'AM172')
# ?mafCompare
# 
# maf1 = read.maf(maf = maf_df[maf_df$Origin_sample_Barcode == 'AM168', ])
# maf2 = read.maf(maf = maf_df[maf_df$Origin_sample_Barcode == 'AM172', ])
# comp_obj <- mafCompare(m1 = maf1, 
#                        m2 = maf2, 
#                        m1Name = 'AM168', m2Name = 'AM172', minMut = 2)
# print(comp_obj)

## AM172 vs AM179

draw_co_oncoplot(maf_df[maf_df$Origin_sample_Barcode == 'AM172', ],
                 maf_df[maf_df$Origin_sample_Barcode == 'AM179', ],
                 'AM172', 'AM179')

## AM168 vs AM179

draw_co_oncoplot(maf_df[maf_df$Origin_sample_Barcode == 'AM168', ],
                 maf_df[maf_df$Origin_sample_Barcode == 'AM179', ],
                 'AM168', 'AM179')

## Fib86 vs Fib172

draw_co_oncoplot(maf_df[maf_df$Origin_sample_Barcode == 'Fib86', ],
                 maf_df[maf_df$Origin_sample_Barcode == 'Fib172', ],
                 'Fib86', 'Fib172')

## Blo172 vs Blo86

draw_co_oncoplot(maf_df[maf_df$Origin_sample_Barcode == 'Blo172', ],
                 maf_df[maf_df$Origin_sample_Barcode == 'Blo86', ],
                 'Blo172', 'Blo86')


# origin group 단위 비교

## AM vs Fib

draw_co_oncoplot(maf_df[maf_df$Origin_class == 'AM', ],
                 maf_df[maf_df$Origin_class == 'Fib', ],
                 'AM', 'Fib')

## AM vs Blo
draw_co_oncoplot(maf_df[maf_df$Origin_class == 'AM', ],
                 maf_df[maf_df$Origin_class == 'Blo', ],
                 'AM', 'Blo')

## Blo vs Fib
draw_co_oncoplot(maf_df[maf_df$Origin_class == 'Blo', ],
                 maf_df[maf_df$Origin_class == 'Fib', ],
                 'Blo', 'Fib')



# Age

## Old vs Fetal
draw_co_oncoplot(maf_df[maf_df$Donor_age == 'Old', ],
                 maf_df[maf_df$Donor_age == 'Fetal', ],
                 'Old', 'Fetal')



# Mt_mut

## mut vs wt
draw_co_oncoplot(maf_df[maf_df$Mt_mutation == 'True', ],
                 maf_df[maf_df$Mt_mutation == 'False', ],
                 'Mt_mut', 'wild_type')

setwd(old.path)
getwd()




































