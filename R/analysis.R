library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)



D <- read_excel("data/data.xlt")
pep_counts = table(D$Peptide)
multiple_peps = names(pep_counts[pep_counts > 1])

# if more than one, how many peptides appeared in more than one case?
table(pep_counts[pep_counts > 1])

# Filtering out multiply occuring peptides.
D2 = D %>% filter(!Peptide %in% multiple_peps)
colnames(D2)[1] = "Prot"
DL = D2 %>% select(Prot, Peptide, contains('-')) %>% 
  gather(des, I, 3:14)

min(DL$I[DL$I > 0])

hist(log10(DL$I[DL$I > 0]), breaks = 1000)
DL = DL %>% separate(des, c('cond', 'repl'), sep='-')

PepCond = DL %>% 
  group_by(Peptide, cond) %>% 
  summarize(n_repl = sum(I > 0),
            med_I = median(I[I > 0]))

ggplot(PepCond, aes(x=ordered(n_repl), y=med_I)) + 
  geom_boxplot() +
  scale_y_log10() + 
  facet_grid(.~cond) +
  xlab('No of replicates peptide was observed in') +
  ylab('median intensity')
# with(PepCond, plot(n_repl, med_i))

# no effect of technical replicate
DL %>% filter(I > 0) %>%
  ggplot(aes(x=repl, y=I)) + geom_boxplot() +
  facet_grid(.~cond) + 
  scale_y_log10()

  
ggplot(PepCond, aes(x=ordered(n_repl), y=med_I)) + 
  geom_boxplot() +
  scale_y_log10() + 
  facet_grid(.~cond) +
  xlab('No of replicates peptide was observed in') +
  ylab('median intensity')



peptides_per_proteins =
  DL %>% 
  group_by(Prot) %>% 
  summarize(n_peptides = n_distinct(Peptide)) %>%
  group_by(n_peptides) %>%
  count %>% ungroup %>% rename(n_proteins=n)

nnmedians =
  DL %>% 
  group_by(Prot, Peptide, cond) %>%
  summarize(
    nn_I_n  = sum(I > 0),
    I_nnmed = median(I[I >0])
  ) %>%
  ungroup %>%
  mutate(I_nnmed = ifelse(is.na(I_nnmed), 0, I_nnmed))

nnmedians_w = 
  nnmedians %>% 
  select(-nn_I_n) %>% 
  spread(cond, I_nnmed)

normalize = function(A, B) ifelse(A+B > 0, (A-B)/(A+B), 0)
nnmedians_w = nnmedians_w %>%
  mutate(
    W.C2_Dark = normalize(C2, Dark),
    W.Light15min_Dark = normalize(Light15min, Dark),
    W.Light2h_Dark = normalize(Light2h, Dark)
  )

all_cond_g1 =
  nnmedians %>% group_by(Prot, Peptide) %>%
  summarize(all_cond_g1 = all(nn_I_n > 0)) %>%
  filter(all_cond_g1) %>%
  select(-all_cond_g1) %>% ungroup

nnmedians_w_cond_g1 = left_join(all_cond_g1, nnmedians_w)

conds = colnames(nnmedians_w)[3:6]
WW = nnmedians_w %>% 
  select(-conds) %>%
  gather(cond, W, 3:5)

WW %>% filter(W != 0, W != -1, W != 1) %>%
  ggplot(aes(x=W)) + 
  geom_histogram(bins=1000)

WW_cont = WW %>% filter(W != 0, W != -1, W != 1)
WW_disc = WW %>% filter(W == 0| W == -1 | W == 1)

WW_cont %>%
  ggplot(aes(x=W)) +
  geom_density() +
  facet_wrap(~cond)

WW_disc %>%
  ggplot(aes(x=W)) + geom_bar() +
  facet_wrap(~cond)
# the only advantage of the above, that there is a natural support in [-1, 1]

## simplify the analysis to more reasonable quantities: the fold-changes and the discrete distribution for: one or more conditions missing.


nnmedians_w = nnmedians_w %>%
  mutate(
    fc.C2_Dark = C2/Dark,
    fc.Light15min_Dark = Light15min/Dark,
    fc.Light2h_Dark = Light2h/Dark
  )

DD = nnmedians_w %>% select(-contains('W.'))
zero_dark = DD %>% filter(Dark == 0)
nonzero_dark = DD %>% filter(Dark > 0)

nonzero_dark_C2 = nonzero_dark %>% filter(C2 > 0)
nonzero_dark_zero_C2 = nonzero_dark %>% filter(C2 == 0)
nonzero_dark_L15min = nonzero_dark %>% filter(Light15min > 0)
nonzero_dark_zero_L15min = nonzero_dark %>% filter(Light15min == 0)
nonzero_dark_L2h = nonzero_dark %>% filter(Light2h > 0)
nonzero_dark_zero_L2h = nonzero_dark %>% filter(Light2h == 0)

zero_dark_nonzero_C2 = zero_dark %>% filter(C2 > 0)
zero_dark_nonzero_L15min = zero_dark %>% filter(Light15min > 0)
zero_dark_nonzero_L2h = zero_dark %>% filter(Light2h > 0)

# X = nonzero_dark_C2
find_significant = function(X){
  XL = X %>% 
    select(-conds) %>%
    gather(cond, fc, 3:5) %>%
    mutate(cond = str_remove(cond, 'fc.'))
  
  XL = XL %>% mutate(log2fc = log2(fc)) %>%
    group_by(cond) %>%
    mutate(logfc_stdev = (log2fc - median(log2fc))/mad(log2fc) ) %>%
    ungroup 
  
  # count how many times a protein was spotted with a fc larger than exp(3)
  UP = XL %>% filter(logfc_stdev > 3) %>% count(Prot)
  DOWN = XL %>% filter(logfc_stdev < -3) %>% count(Prot)
  list(XL=XL, UP=UP, DOWN=DOWN)
}

nonzero_dark_C2_sig = find_significant(nonzero_dark_C2)
nonzero_dark_L15min_sig = find_significant(nonzero_dark_L15min)
nonzero_dark_L2h_sig = find_significant(nonzero_dark_L2h)

nonzero_dark_C2_sig$XL %>%
  ggplot(aes(x=fc, color=cond)) + geom_density() + scale_x_log10()

# what's the intersection?
intersect(unique(nonzero_dark_C2_sig$UP$Prot), unique(nonzero_dark_C2_sig$DOWN$Prot))
# what to do with those?

# what to do with 0s now? Assume, they are significant.
# possible problem: we treat all these the same, as upregulated
nonzero_dark_zero_C2 %>% filter(Prot == "TGGT1_201760-t26_1-p1")

nonzero_dark = bind_rows(
    nonzero_dark_zero_C2 %>% count(Prot) %>% mutate(cond = 'C2'),
    nonzero_dark_zero_L15min %>% count(Prot) %>% mutate(cond = 'Light15min'),
    nonzero_dark_zero_L2h %>% count(Prot) %>% mutate(cond = 'Light2h')
  ) %>% 
  mutate(reg = "DOWN")

zerodark = bind_rows(
    zero_dark_nonzero_C2 %>% count(Prot) %>% mutate(cond = 'C2'),
    zero_dark_nonzero_L15min %>% count(Prot) %>% mutate(cond = 'Light15min'),
    zero_dark_nonzero_L2h %>% count(Prot) %>% mutate(cond = 'Light2h')
  ) %>% 
  mutate(reg = "UP")

# neglect_zero_zero??? not observed anywhere, but some other condition.
zero_dark_or_other = bind_rows(nonzero_dark, zerodark)
zero_dark_or_other$fc = FALSE


nonzeros = bind_rows(
  nonzero_dark_C2_sig$UP %>% mutate(cond='C2', reg='UP', fc=T),
  nonzero_dark_L15min_sig$UP %>% mutate(cond='Light15min', reg='UP', fc=T),
  nonzero_dark_L2h_sig$UP %>% mutate(cond='Light2h', reg='UP', fc=T),
  nonzero_dark_C2_sig$DOWN %>% mutate(cond='C2', reg='DOWN', fc=T),
  nonzero_dark_L15min_sig$DOWN %>% mutate(cond='Light15min', reg='DOWN', fc=T),
  nonzero_dark_L2h_sig$DOWN %>% mutate(cond='Light2h', reg='DOWN', fc=T)
)
nonzeros %>% count(Prot)
nonzeros %>%
  select(-fc) %>%
  spread(reg, n, fill=0) %>%
  count(DOWN, UP) %>%
  ggplot(aes(x=DOWN, y=UP, size=n)) + geom_point()

ultimate = bind_rows(zero_dark_or_other, nonzeros)
U = ultimate %>% group_by(Prot, reg) %>% summarize(n = sum(n)) %>% ungroup
U %>% spread(reg, n, fill=0) %>%
  count(DOWN, UP) %>%
  ggplot(aes(x=DOWN, y=UP, size=n)) + geom_point()

U %>% filter(n==50)
D %>% filter(D[,1] == 'TGGT1_311230-t26_1-p1') %>% 

