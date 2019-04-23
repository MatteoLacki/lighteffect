library(readxl)
pvals <- read_excel("data/MASTERFILE-worked on it - TM-OG-V1.xlsx", 
                    sheet = "p values")

pvals = pvals %>% select(contains('ttest'))

methods = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
R = lapply(
  methods,
  function(m) apply(pvals, 2, p.adjust, method=m) 
)
names(R) = methods

res = sapply(methods,
       function(m, alpha=.05) colSums(R[[m]] < alpha/4, na.rm=T))

write.csv(R$fdr, file="data/fdr.csv")
write.csv(R$fdr * 4, file="data/fdr4.csv")

