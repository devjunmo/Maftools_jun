
getwd()

library(dplyr)
library(tibble)

name <- c('kim','lee', 'lee', 'park', 'choi')
age <- c(20, 18, 22, 30, 40)
grade <- c('H','H', 'L', 'L', 'L')
Group <- c('a', 'a','a', 'a', 'a')

name2 <- c('kim', 'lee', 'lee', 'park', 'seo')
age2 <- c(20,21,19, 18, 30)
grade2 <- c('H','H','H', 'H', 'L')
Group2 <- c('b','b', 'b','b', 'b')

df1 <- tibble(data.frame(name, age, grade, Group))

df2 <- tibble(data.frame(name2, age2, grade2, Group2))

names(df2) <- c('name','age', 'grade', 'Group')

df1
df2

# full.join.tst <- full_join(df1, df2, by=c('name', 'age'))
full.join.tst <- full_join(df1, df2, by=c('name', 'age', 'grade'))
# full.join.tst <- full_join(df1, df2, by=c('name', 'age', 'grade', 'Group')) # 공통인건 하나만, 다른것들은 모두 나열
# left.join.tst <- left_join(df1, df2, by=c('name', 'age', 'grade', 'Group'))
                                                                            # 이새끼랑 df1을 차집합해서 이새끼 sp을 구하면 df2 sp가 뜨겠네
# full.join.tst <- full_join(df1, df2, by=c('name', 'grade'))
# full.join.tst <- full_join(df1, df2, by=c('name'))
full.join.tst
left.join.tst

full.join.tst[full.join.tst$name=='lee', ]

print(sum(is.na(full.join.tst))) # 결측치 존재 확인
print(colSums(is.na(full.join.tst))) # 컬럼별로 NA 확인

inc.na <- full.join.tst[as.logical(apply(is.na(full.join.tst), MARGIN = 1, FUN = sum)), ]

inc.na$name
























