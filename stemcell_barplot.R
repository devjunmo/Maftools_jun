getwd()

setwd('C:/Users/user/Desktop/git_manager/Maftools_jun')
getwd()


library(tibble)
library(dplyr)

root_dir = 'E:/stemcell/somatic_analysis/maf/mutect2/DP_AF_filtered_maf/exonic_maf'
old.path <- setwd(root_dir)

# clinl.info = './Tera_comp/Tera_cohort_info_220321.csv'
clinl.info = './cohort_comp/IPS_cohort_info_220321.csv' # ips, heso

getwd()

# maf file list 가져오기

f_lst <- list.files(pattern = '*.maf')

mutation_count_lst <- list()
mutation_count_lst
class(mutation_count_lst)

# mutation_count_lst[['test']] <- 651
# mutation_count_lst


# mutation count df 만들기

samplesVec <- c()
mutCount <- c()

for(file in f_lst){
  
  f_name <- strsplit(file, '_')[[1]][1]
  print(f_name)
  samplesVec <- c(samplesVec, f_name)
  mut_count <- dim(read.delim(file))[1]
  mutCount <- c(mutCount, mut_count)
  
}

samplesVec
mutCount


mutation_df <- tibble(data.frame(samplesVec, mutCount))

mutation_df <- rename(mutation_df, 
                      Tumor_Sample_Barcode=samplesVec,
                      Mutation_Count=mutCount)

mutation_df

# clin.info 가져오기

clin.info.df <- tibble(read.csv(clinl.info, fileEncoding="UTF-8-BOM",
                                na.strings = c('', 'NA')))

print(clin.info.df, width = Inf) # 28 x 6

clin.info.df <- select(clin.info.df, 
                       Tumor_Sample_Barcode, Origin_sample, 
                       Origin_class, Donor_age, Mt_mutation)

print(clin.info.df, width = Inf, n=Inf) # 28 x 6


# NA 제거 

if(sum(is.na(clin.info.df$Tumor_Sample_Barcode))!=0){
  clin.info.df <- clin.info.df[!is.na(clin.info.df$Tumor_Sample_Barcode), ]
}


print(clin.info.df, width = Inf, n=Inf)


# maf, clin.info inner join

mutation_df # 61 x 2
clin.info.df # 22 x 6


joined_maf <- inner_join(clin.info.df, mutation_df,
                         by=c('Tumor_Sample_Barcode'))

print(joined_maf, n =Inf)

unique(joined_maf$Origin_sample)


maf_class <- unique(joined_maf$Origin_class)
maf_class

joined_maf_1 <- arrange(joined_maf,
                      Origin_class, Origin_sample)
joined_maf_1

joined_maf_1 <- mutate(joined_maf_1, Origin_class_plot=Origin_class)
print(joined_maf_1, n =Inf)

joined_maf_1$Origin_class <- as.factor(joined_maf_1$Origin_class)
print(joined_maf_1, n =Inf)


joined_maf_1 <- transform(joined_maf_1,
                          Origin_class=factor(Origin_class, 
                                               levels = unique(joined_maf_1$Origin_class)))


joined_maf_1$Origin_class

levels(joined_maf_1$Origin_class) <- rep(1:length(unique(joined_maf_1$Origin_class)))
joined_maf_1$Origin_class


library(ggplot2)
joined_maf_1$Origin_class <- as.numeric(joined_maf_1$Origin_class)

paste(joined_maf_1$Mutation_Count, joined_maf_1$Origin_sample,
      sep = '\n')

joined_maf_1$annotTxt <- paste(joined_maf_1$Mutation_Count, joined_maf_1$Origin_sample,
                               sep = '\n')


ggplot(joined_maf_1, aes(x=reorder(Tumor_Sample_Barcode, Origin_class),
                       y=log(Mutation_Count),
                       fill=Origin_class_plot)) +
  geom_bar(stat="identity", position = 'dodge') + 
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  geom_text(aes(label=c(annotTxt)), vjust=1.6, color="black", size=3.5)




# 
# 
# joined_maf_2 <- arrange(joined_maf,
#                         Origin_class, Origin_sample)
# 
# joined_maf_2 <- mutate(joined_maf_2, Origin_samples_plot=Origin_sample)
# 
# print(joined_maf_2, n =Inf)
# 
# joined_maf_2$Origin_sample <- as.factor(joined_maf_2$Origin_sample)
# 
# print(joined_maf_2, n =Inf)
# 
# # joined_maf_2$Origin_sample <- rep(1:length(unique(joined_maf_2$Origin_sample)))
# joined_maf_2$Origin_sample
# levels(joined_maf_2$Origin_sample) <- rep(1:length(unique(joined_maf_2$Origin_sample)))
# print(joined_maf_2, n =Inf)
# sort(joined_maf_2$Origin_sample)
# joined_maf_2$Origin_sample <- sort(joined_maf_2$Origin_sample)
# print(joined_maf_2, n =Inf)
# 
# # joined_maf_2$Origin_sample <- sort(joined_maf_2$Origin_sample)
# # print(joined_maf_2, n =Inf)
# 
# library(ggplot2)
# joined_maf_2$Origin_sample <- as.numeric(joined_maf_2$Origin_sample)
# # joined_maf_2$Origin_samples_plot <- as.character(joined_maf_2$Origin_sample)
# print(joined_maf_2, n =Inf)
# 
# 
# ggplot(joined_maf_2, aes(x=reorder(Tumor_Sample_Barcode, Origin_sample),
#                          y=log(Mutation_Count),
#                          fill=Origin_samples_plot)) +
#   geom_bar(stat="identity", position = 'dodge') +
#   theme(axis.text.x=element_text(angle=60, hjust=1)) +
#   geom_text(aes(label=Mutation_Count), vjust=1.6, color="black", size=3.5) +
#   scale_color_discrete(breaks=unique(joined_maf_2$Origin_samples_plot))






# joined_maf_2 <- arrange(joined_maf,
#                         Origin_class, Origin_sample)
# 
# 
# 
# joined_maf_2 <- mutate(joined_maf_2, Origin_samples_plot=Origin_sample)
# 
# print(joined_maf_2, n =Inf)
# 
# joined_maf_2$Origin_sample <- as.factor(joined_maf_2$Origin_sample)
# 
# print(joined_maf_2, n =Inf)
# 
# joined_maf_2 <- transform(joined_maf_2,
#                           Origin_sample=factor(Origin_sample, 
#                                                levels = unique(joined_maf_2$Origin_sample)))
# joined_maf_2
# joined_maf_2$Origin_sample
# 
# levels(joined_maf_2$Origin_sample) <- rep(1:length(unique(joined_maf_2$Origin_sample)))
# joined_maf_2$Origin_sample
# joined_maf_2
# 
# 
# # joined_maf_2$Origin_sample <- sort(joined_maf_2$Origin_sample)
# # print(joined_maf_2, n =Inf)
# 
# library(ggplot2)
# joined_maf_2$Origin_sample <- as.numeric(joined_maf_2$Origin_sample)
# # joined_maf_2$Origin_samples_plot <- as.character(joined_maf_2$Origin_sample)
# # print(joined_maf_2, n =Inf)
# 
# 
# ggplot(joined_maf_2, aes(x=reorder(Tumor_Sample_Barcode, Origin_sample),
#                          y=log(Mutation_Count),
#                          fill=Origin_samples_plot)) +
#   geom_bar(stat="identity", position = 'dodge') +
#   theme(axis.text.x=element_text(angle=60, hjust=1)) +
#   geom_text(aes(label=Mutation_Count), vjust=1.6, color="black", size=3.5)


# ?scale_color_discrete
# 
# label.tst <- paste0("Tree",1:5)
# class(label.tst)
# class(unique(joined_maf_2$Origin_samples_plot))
# 





# age에 따라 mutation plotting

print(joined_maf, n=Inf)

joined_maf_age <- arrange(joined_maf, Donor_age)
print(joined_maf_age, n=Inf)

joined_maf_age <- joined_maf_age[!joined_maf_age$Donor_age=='amb',]
print(joined_maf_age, n=Inf)



# joined_maf_2$Origin_sample <- sort(joined_maf_2$Origin_sample)
# print(joined_maf_2, n =Inf)

library(ggplot2)


ggplot(joined_maf_age, aes(x=Tumor_Sample_Barcode,
                         y=log(Mutation_Count),
                         fill=Donor_age)) +
  geom_bar(stat="identity", position = 'dodge') +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  geom_text(aes(label=Mutation_Count), vjust=1.6, color="black", size=3.5)





# mt mutations 에 따라 mutation plotting


print(joined_maf, n=Inf)

joined_maf_mtDNA <- arrange(joined_maf, Mt_mutation)
print(joined_maf_mtDNA, n=Inf)
na.omit(joined_maf_mtDNA)

joined_maf_mtDNA_na.omit <- na.omit(joined_maf_mtDNA)
print(joined_maf_mtDNA_na.omit, n=Inf)
joined_maf_mtDNA_na.omit <- mutate(joined_maf_mtDNA_na.omit,
                                   mt_for_plot=Mt_mutation)

joined_maf_mtDNA_na.omit$mt_for_plot <- 
  as.numeric(joined_maf_mtDNA_na.omit$mt_for_plot)

print(joined_maf_mtDNA_na.omit, n=Inf)


ggplot(joined_maf_mtDNA_na.omit, aes(x=Tumor_Sample_Barcode,
                                     y=log(Mutation_Count),
                                     fill=Mt_mutation)) +
  geom_bar(stat="identity", position = 'dodge') +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  geom_text(aes(label=Mutation_Count), vjust=1.6, color="black", size=3.5)


ggplot(joined_maf_mtDNA_na.omit, aes(x=reorder(Tumor_Sample_Barcode, mt_for_plot),
                           y=log(Mutation_Count),
                           fill=Mt_mutation)) +
  geom_bar(stat="identity", position = 'dodge') +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  geom_text(aes(label=Mutation_Count), vjust=1.6, color="black", size=3.5)





joined_maf_mtDNA_mean_exceptCancer <- joined_maf_mtDNA_na.omit %>%
  filter(Origin_class!='Cancer') %>%
  group_by(Mt_mutation) %>%
  summarise(mean_mutations = mean(Mutation_Count))

joined_maf_mtDNA_mean_exceptCancer

print(joined_maf_mtDNA_na.omit.filter, n=Inf)


ggplot(joined_maf_mtDNA_mean_exceptCancer, 
       aes(x=Mt_mutation, y=mean_mutations, fill=Mt_mutation)) +
  geom_bar(stat="identity", position = 'dodge') +
  geom_text(aes(label=round(mean_mutations, digits = 2)), vjust=1.6, color="black", size=3.5)
























