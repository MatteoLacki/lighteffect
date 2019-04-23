library(readr)
library(gplots)
PPP <- read_delim("data/ClustVis_Heatmap_V1.csv", 
                  ";", 
                  escape_double = FALSE, 
                  locale = locale(decimal_mark = ",", 
                                  grouping_mark = "."),
                  trim_ws = TRUE)

prot2pep <- read_excel("data/Protein-Peptide Name ALL.xlsx")
colnames(prot2pep)[1] = 'Protein'
PPP = left_join(PPP, prot2pep)
PPPM = data.matrix(PPP[,4:7])
my_palette <- colorRampPalette(c("red", "black", "green"))(n = 40)

rownames(PPPM) = PPP$Peptide
pdf('plots/hm0.pdf',    # create PNG for the heat map        
    width = 10,        # 5 x 300 pixels
    height = 30)        # smaller font size
heatmap.2(PPPM, 
          lhei=c(1,10),
          trace = 'none',
          col=my_palette,
          margins =c(12,20))
dev.off()       


rownames(PPPM) = PPP$Protein
pdf('plots/hm1.pdf',    # create PNG for the heat map        
    width = 10,        # 5 x 300 pixels
    height = 30)        # smaller font size
heatmap.2(PPPM, 
          lhei=c(1,10),
          trace = 'none',
          col=my_palette,
          margins =c(12,20))
dev.off()       

rownames(PPPM) = PPP$Peptide
PPPM = PPPM[,-4]
pdf('plots/hm0_nolight.pdf',    # create PNG for the heat map        
    width = 10,        # 5 x 300 pixels
    height = 30)        # smaller font size
heatmap.2(PPPM, 
          lhei=c(1,10),
          trace = 'none',
          col=my_palette,
          margins =c(12,20))
dev.off()       


rownames(PPPM) = PPP$Protein
pdf('plots/hm1_nolight.pdf',    # create PNG for the heat map        
    width = 10,        # 5 x 300 pixels
    height = 30)        # smaller font size
heatmap.2(PPPM, 
          lhei=c(1,10),
          trace = 'none',
          col=my_palette,
          margins =c(12,20))
dev.off()       
