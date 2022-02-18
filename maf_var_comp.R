
# maf 파일간 sample id를 비교하여 동일한 sample id에 동일한 변이가 있는지 확인

old.path <- setwd('C:/Users/user/Desktop/Asannet_data/wd')
getwd()

library(dplyr)
library(tibble)
library(ggplot2)

old.maf <- tibble(read.delim('panc_processed.maf'))
old.maf$grp <- 'old'
old.maf # 5,647 x 134

length(unique(old.maf$Tumor_Sample_Barcode)) # 687

new.maf.K <- tibble(read.csv('Pancreatic_cancer_K_new.csv'))
new.maf.J <- tibble(read.csv('Pancreatic_cancer_J_new.csv'))

new.maf.K$grp <- 'new.k'
new.maf.J$grp <- 'new.J'



new.maf.K # 5,589 x 38
length(unique(new.maf.K$Tumor_Sample_Barcode)) # 679


new.maf.J # 1,533 x 38
length(unique(new.maf.J$Tumor_Sample_Barcode)) # 190

colnames(new.maf.K)


def_select.columns <- function(maf_, target_col_){
  
  maf_ <- maf_ %>% select(target_col_)
  
  return(maf_)
}



select.col.name <- c(
  "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
  "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "grp"
)

old.maf <- def_select.columns(old.maf, select.col.name)
new.maf.K <- def_select.columns(new.maf.K, select.col.name)
new.maf.J <- def_select.columns(new.maf.J, select.col.name)

def_comp.sam.id <- function(maf1_, maf2_){
  
  maf1.id <- maf1_$Tumor_Sample_Barcode
  maf2.id <- maf2_$Tumor_Sample_Barcode
  
  print("maf1 unique TumorID count")
  print(length(unique(maf1.id)))
  
  print("maf2 unique TumorID count")
  print(length(unique(maf2.id)))
  
  print("Intersect length")
  print(length(intersect(maf1.id, maf2.id)))
  
  return(intersect(maf1.id, maf2.id))
  
}


K.ids <- def_comp.sam.id(old.maf, new.maf.K) # 교집합 모두 만족
J.ids <- def_comp.sam.id(old.maf, new.maf.J) # 교집합 모두 만족 190개



def_full.outer.join <- function(maf1_, maf2_, by_){
  
  full.join.maf <- NA
  full.join.maf <- full_join(maf1_, maf2_, by=by_)
  print(full.join.maf)
  print('NA count')
  print(sum(is.na(full.join.maf))) # 결측치 존재 확인
  print(colSums(is.na(full.join.maf))) # 컬럼별로 NA 확인
  
  return(full.join.maf)
  
}

old.maf
new.maf.K

by=c('Hugo_Symbol', 'Chromosome', 'Start_Position', 
     'End_Position', 'Reference_Allele', 
     'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode')

K.fj.maf <- def_full.outer.join(old.maf, new.maf.K, by=by)
J.fj.maf <- def_full.outer.join(old.maf, new.maf.J, by=by)


def_get_common <- function(maf_){
  
  common_part_maf <- filter(maf_, !is.na(grp.x) & 
                              !is.na(grp.y)) # 샘플단위로 변이가 모두일치
  print(common_part_maf)
  return(common_part_maf)
}

def_get_target_specific <- function(maf_){
  
  sp_part_maf <- filter(maf_, is.na(grp.x))
  return(sp_part_maf)
  
}

def_get_old_specific <- function(maf_){
  
  sp_part_maf <- filter(maf_, is.na(grp.y))
  return(sp_part_maf)
  
}


def_show_specific_var_per_samples <- function(outer_join_maf_){
  
  print('---------------------------------------------------------------------')
  # common
  print('Common part variant maf')
  common_var_maf <- def_get_common(outer_join_maf_)
  match_all.sample <- unique(common_var_maf$Tumor_Sample_Barcode)
  print('common var maf의 unique samples')
  print((match_all.sample))
  print(length(match_all.sample)) # 679
  print('---------------------------------------------------------------------')
  
  print('---------------------------------------------------------------------')
  # specific 1: Old specific
  print('Old part variant maf')
  old.sp.var.maf <- def_get_old_specific(outer_join_maf_)
  
  old.sp.var.maf.samples <- unique(old.sp.var.maf$Tumor_Sample_Barcode)
  print(old.sp.var.maf.samples)
  print(length(old.sp.var.maf.samples))
  print('---------------------------------------------------------------------')
  
  print('---------------------------------------------------------------------')
  # specific 2: target specific
  print('Target part variant maf')
  target.sp.var.maf <- def_get_target_specific(outer_join_maf_)
  target.sp.var.maf.samples <- unique(target.sp.var.maf$Tumor_Sample_Barcode)
  print(target.sp.var.maf.samples)
  print(length(target.sp.var.maf.samples))
  print('---------------------------------------------------------------------')
}


def_show_specific_var_per_samples(K.fj.maf)
def_show_specific_var_per_samples(J.fj.maf)

# K.fj.maf[is.na(K.fj.maf),]


k.fj.nonNA <- filter(K.fj.maf, !is.na(Hugo_Symbol.y)) # 샘플단위로 변이가 모두일치
k.fj.nonNA.match_all.sample <- unique(k.fj.nonNA$Tumor_Sample_Barcode)
k.fj.nonNA.match_all.sample
length(k.fj.nonNA.match_all.sample) # 679

k.fj.NA <- filter(K.fj.maf, is.na(Hugo_Symbol.y))
k.fj.NA.match_all.sample <- unique(k.fj.NA$Tumor_Sample_Barcode)
length(k.fj.NA.match_all.sample)
k.fj.NA.match_all.sample


intersect(unique(k.fj.nonNA$Tumor_Sample_Barcode),
          unique(k.fj.NA$Tumor_Sample_Barcode)) # 겹치는거 X





j.fj.nonNA <- filter(J.fj.maf, !is.na(Hugo_Symbol.y))
j.fj.nonNA.match_all.sample <- unique(j.fj.nonNA$Tumor_Sample_Barcode)
j.fj.nonNA.match_all.sample
length(j.fj.nonNA.match_all.sample) # 679



j.fj.NA <- filter(J.fj.maf, is.na(Hugo_Symbol.y))

intersect(unique(j.fj.nonNA$Tumor_Sample_Barcode),
          unique(j.fj.NA$Tumor_Sample_Barcode)) # 겹치는 거 X




# maf.K.NA <- filter(full.join.test, is.na(Hugo_Symbol.y)) # 조건에 맞는 행출력
# maf.K.NA
# is.na(J.fj.maf)
# unique(maf.K.NA$Tumor_Sample_Barcode)
# print(maf.K.NA, n=nrow(maf.K.NA))
# unique(maf.K.NA$Hugo_Symbol.x)


def_get_y_specific <- function(full.join.df_){
  
  y.NA.maf <- filter(full.join.df_, is.na(Hugo_Symbol.y)) # 조건에 맞는 행출력
  print('---------------------------------------------------------------------')
  print('| y.NA.maf |')
  print('---------------------------------------------------------------------')
  print(y.NA.maf)
  print('---------------------------------------------------------------------')
  print('| y.NA.maf - unique Sample ID |')
  print(unique(y.NA.maf$Tumor_Sample_Barcode))
  print(length(unique(y.NA.maf$Tumor_Sample_Barcode)))
  print('---------------------------------------------------------------------')
  # print(maf.K.NA, n=nrow(maf.K.NA))
  print('| y.NA.maf - unique Gene ID |')
  print(unique(y.NA.maf$Hugo_Symbol.x))
  print('---------------------------------------------------------------------')
  
  return(y.NA.maf)
  
}


K.y.NA.maf <- def_get_y_specific(K.fj.maf) 
J.y.NA.maf <- def_get_y_specific

K.y.NA.maf # old specific한거


# apply(is.na(full.join.test), MARGIN = 1, FUN = sum)!=0
# 
# inc.na <- full.join[apply(is.na(full.join.test), MARGIN = 1, FUN = sum)!=0, ]

# 그림: unqiue로 gene list 뽑고, count로 조지기
?geom_histogram

# 김성주샘
ggplot(K.y.NA.maf, aes(x=Hugo_Symbol.x)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = 'old vs new maf file - 김성주',
       x='Genes - only in old maf files',
       y='Gene count')


# j.hugo.table <- table(J.y.NA.maf$Hugo_Symbol.x)
# class(table(J.y.NA.maf$Hugo_Symbol.x)>1)
# length(j.hugo.table)
# j.hugo.table[301]
# j.hugo.table[[301]]
# class(j.hugo.table[1])
# 
# 
# 
# as.vector(table(J.y.NA.maf$Hugo_Symbol.x))
# table(J.y.NA.maf$Hugo_Symbol.x)>1
# as.vector(table(J.y.NA.maf$Hugo_Symbol.x)>1)


# J.y.NA.maf[as.vector(table(J.y.NA.maf$Hugo_Symbol.x)>1)]


J.y.NA.maf # 4,113 x 13
J.y.NA.maf[duplicated(J.y.NA.maf$Hugo_Symbol.x),]
dim(J.y.NA.maf[duplicated(J.y.NA.maf$Hugo_Symbol.x),]) # 3812   13

J.y.cut.onlyone_gene <- J.y.NA.maf[duplicated(J.y.NA.maf$Hugo_Symbol.x),]
J.y.cut.onlyone_gene[J.y.cut.onlyone_gene$Hugo_Symbol.x=='APC',]
?theme

ggplot(J.y.cut.onlyone_gene, aes(x=Hugo_Symbol.x)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6, face='bold')) +
  labs(title = 'old vs new maf file - 조의리',
       x='Genes - only in old maf files',
       y='Gene count')



print('end')
setwd(old.path)
getwd()











