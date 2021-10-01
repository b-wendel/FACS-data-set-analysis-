library(tidyverse)
library(ggplot2)
library(ggpubr)
library(plotly)
library(ggforce)
library(readr)
library(rebus)
library(stringr)
library(hrbrthemes)
library(ggtext)
library(rstatix)
library(gridExtra)
library(viridis)
library(grid)
#install.packages("pheatmap")
library(pheatmap)
#install.packages("GGally")
library(GGally)
library(RColorBrewer)
library(writexl)
library(Hmisc)
library(reshape2)
library(fdrtool)

library(corrplot)



setwd("~/Ped-Covid Studie")
load("~/Ped-Covid Studie/FACS/scripts/clinical_table_ws_full.RData")
load("~/Ped-Covid Studie/FACS/scripts/all_panels_common_WS.RData")
#load("~/Ped-Covid Studie/FACS/scripts/correlation_all_panels_clinical_WS_full.RData")

complete_panels <- Reduce(intersect, list(panel1_cov$merge_ID,panel2_cov$merge_ID,P3_cov_complete$merge_ID))

p1_common <- panel1_cov[panel1_cov$merge_ID %in% complete_panels,]
p2_common <- panel2_cov[panel2_cov$merge_ID %in% complete_panels,]
p3_common <- P3_cov_complete[P3_cov_complete$merge_ID %in% complete_panels,]

#merge all panels together, exclude outliers , make second timepoint of COV068 day 0 because PCR tst was positive for the first time only before that
pall <- p1_common[,c(6:10,12,13,15,17:19,21,22,25:41,52,53,55,57,58,60,61,63,64,66,67,76,77,79,80,172)] %>% 
  left_join(p2_common[,c(6:15,98)], by = "merge_ID") %>% 
  left_join(p3_common[,c(4:19,24:29,33,116)], by = "merge_ID") %>% 
  left_join(p1_common[,90:181], by = "merge_ID")



pall <- pall %>% 
  mutate(Disease_course = cut(max_WHO_classification, breaks=c(-Inf, 0.9, 1.9, 3.9, 5.9,Inf), 
                                             label=c("Healthy", "Asymptomatic", "Mild", "Moderate","Severe"))) %>% 
  mutate(Disease_course = factor(Disease_course, levels = c("Healthy","Asymptomatic","Mild","Moderate","Severe"))) %>% 
  mutate(Age_group = cut(Age, breaks=c(-Inf, 0.9, 10.9, 18.9,Inf), 
                         label=c("<1", "1-10", "11-18", ">18"))) %>% 
  mutate(Age_group = factor(Age_group, levels = c("<1","1-10","11-18",">18"))) %>% 
  mutate(Status = factor(Status , levels=c("Healthy", "COVID", "Non-COVID","MISC"))) %>%
  mutate(days_symptom_start = (pall$sampling_date - pall$date_of_1st_symptoms)) %>% 
  mutate(days_symptom_start = as.numeric(str_remove(days_symptom_start, " days"))) %>% 
  mutate(sex01 =  pall$`Sex_(f=0;_m=1)`) %>% 
  mutate(COVID = pall$covid <- ifelse(pall$Status == "COVID", 1, 0))

#generate combined dataframe 
pall_corr <- pall %>% filter(outlier_code <5 & Age > 0 &  day ==0 & 
                               Status != "Non-COVID")


lab <- read_delim("~/Ped-Covid Studie/metadaten/20210517_lab.csv", 
                  ";", escape_double = FALSE, col_types = cols(.default = "n", 
                                                               analysis.identifier = col_character(),
                                                               Comments= col_character()), 
                  locale = locale(decimal_mark = ","), 
                  trim_ws = TRUE)
#only relevant columns 
colnames(lab)
lab_clean <- lab[,c(1,2,5,8,19,22,26,30,42,45,52,59,62,65,68,71,74,77,83,86,89,95,98,104,110,113,116,119, 142,144)]
colnames(lab_clean) <- str_replace_all(colnames(lab_clean), coll(" (on admission)"), "")

lab_clean <- lab_clean %>% left_join(just_meta[,c(1:25)], by = "analysis.identifier")
lab_clean$Status <- factor(lab_clean$Status , levels=c("Healthy", "Non-COVID", "COVID", "MISC"))

lab_clean <- lab_clean %>% filter(is.na(Status) ==FALSE)



lab_corr <- as.data.frame(lab_clean[lab_clean$analysis.identifier %in% pall_corr$analysis.identifier.x,])
pall_corr <- as.data.frame(pall_corr[pall_corr$analysis.identifier.x %in% lab_corr$analysis.identifier,]) 

lab_corr$analysis.identifier.x <- lab_corr$analysis.identifier


combined <- merge(pall_corr, lab_corr, by = "analysis.identifier.x")
rownames(combined) <- combined$analysis.identifier.x
colnames(combined)

#combined <- combined %>% filter(analysis.identifier.x != "COV006")

correlation_canon <- c("CD3","CD4",                       "CD8" ,                     
                       "Th2",                       "Th17",                      "Th1_CCR4neg",              
                       "CD38_CD4m" , "activ_CD4m",   "Foll_Th_like"    ,         
                       "Tregs",         "EM_CD4"    ,   "CM_CD4"     ,               "naive_CD4" ,
                       "EM_CD8",                    "CM_CD8"  ,                  "naive_CD8" ,
                       "TEM_CD8" ,          "CD38_CD8m" ,                "activ_CD8m"  ,          
                       "CD19","CD27_IgD_DN",                      "transitional",              "plasmablasts"     ,        
                       "naive_CD19"           ,     "non_switch_CD19m"  ,        "class_switch_CD19m"  ,     
                       "CD38_CD21_DN"        ,      "immat_CD19",  "classical_NK"    ,
                       "CD57pos_NK"      ,          "undiff_NK"  ,             
                       "NKTs"   ,                   "ab_CD3"    ,                "yd_CD3",                   
                       "CD4_senesc" ,               "CD8_senesc",                "monocytes",                
                       "non_class_mono"  ,          "intermed_mono",             "class_mono",               
                       "mDC",                       "pDC" ,                     "CD14_CD56_DP" ,
                       "COVID", "max_WHO_classification.x","fever",
                       "days_symptom_start"    , "Supplemental_oxygen",                                                    
                       "Age.x","sex01" ,"BMI.x"  , "PEC_(y=1;n=0)"     ,                                                            
                       "Hemoglobin g/dl" , "Thrombozyten G/l" ,         "Leukozyten G/l"      ,
                       "Neutrophile %"            ,
                       "Lymphocytes %"        ,     "Monocytes %"  ,             "CRP mg/dl",                
                       "Procalcitonin ng/ml" ,      "Ferritin mg/dl",            "IL-6 (<5,9 pg/ml)",        
                       "Kreatinin mg/dl"  ,         "Harnstoff mg/dl" ,          "LDH U/l",                  
                       "Quick %"            ,       "pTT sec"      ,             "Fibrinogen mg/dl"  ,       
                       "D-Dimer ym/ml"          ,   "Antithrombin %"   )

combined <- combined[,correlation_canon] %>% 
  mutate(Supplemental_oxygen = as.numeric(Supplemental_oxygen)) %>% 
  mutate(Supplemental_oxygen = replace(Supplemental_oxygen, is.na(Supplemental_oxygen) == TRUE, 0))



# correlate using pairwise complete observations 
result<-rcorr(as.matrix(combined), type="spearman")
result$n %>% view()
#result$r[result$n<7]<-NA # ignore less than 7 observations
rownames(result$r)
result$P %>% view()


#extract relevant columns & rows c
corrmat <- result$r[44:70,1:43]
signmat <- result$P[44:70,1:43]
#fdrtool(melted_signmat$P, statistic = "pvalue",cutoff.method="fndr")

#add statistical significance levels 
melted_cormat <- melt(corrmat)
melted_cormat$merge <- paste(melted_cormat$Var1, melted_cormat$Var2, sep = "_")
melted_signmat <-  melt(signmat)
melted_signmat$merge <- paste(melted_signmat$Var1, melted_signmat$Var2, sep = "_")
colnames(melted_signmat)[3] <- "P"

melted_cormat <- merge(melted_cormat, melted_signmat, by = "merge")
melted_cormat <- melted_cormat[,c(1:4,7)]
melted_cormat$stars <- cut(melted_cormat$P, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) ## cut??

melted_cormat_fdr <- fdrtool(melted_cormat$P, statistic = "pvalue",cutoff.method="fndr")
melted_cormat_fdr$pval
length(melted_cormat_fdr$lfdr[melted_cormat_fdr$lfdr<.1])

melted_cormat$lfdr <-melted_cormat_fdr$lfdr 

fdrmat <- reshape(melted_cormat[,c(2,3,7)], idvar = "Var1.x", timevar = "Var2.x", direction= "wide")
rownames(fdrmat) <- fdrmat$Var1.x
fdrmat <- fdrmat[,-1]
colnames(fdrmat) <- str_remove(colnames(fdrmat), "lfdr.")
fdrmat <- fdrmat[rownames(corrmat),colnames(corrmat)]

rownames(corrmat) <- c("COVID",                    "max. WHO classification", "Fever",                   
                       "Days since 1st Symptom"  ,     "Supplemental oxygen"  ,    "Age",                   
                       "Male Sex"       ,             "BMI"  ,                  "Comorbidities" ,          
                       "Hemoglobin",          "Platelets" ,        "Leukocytes",          
                       "Neutrophils",            "Lymphocytes",            "Monocytes",             
                       "CRP",                "Procalcitonin",      "Ferritin",          
                       "IL-6",        "Creatinine",          "Urea",         
                       "LDH",                  "Quick",                  "pTT",                 
                       "Fibrinogen",         "D-Dimer",            "Antithrombin") 
colnames(corrmat) <- c("T cells",                "CD4 T cells"   ,             "CD8 T cells" ,               "Th2"  ,             
                       "Th17",               "Th1",        "CD38+ CD4",          "Active CD4",        
                       "fTh",       "Tregs"  ,            "EM CD4 "           ,  "CM CD4"  ,          
                       "Naive CD4",          "EM CD8",             "CM CD8",              "Naive CD8",         
                       "TEM CD8",             "CD38+ CD8" ,         "Active CD8" ,        "B cells"  ,            
                       "CD27 IgD DN B cells",        "Transitional B cells",       "Antibody-secreting cells",       "Naive B cells",        
                       "Non-switch B cells",   "Class-switch B cells" ,"CD38 CD21 DN B cells" ,      "Immatature B cells",        
                       "Classical NK",       "CD57+ NK",         "Undifferentiated NK",          "NKT",              
                       "ab-TCR T cells"   ,          "yd-TCR T cells",             "Senescent CD4",        "Senescent CD8",        
                       "Monocytes"    ,      "Non-classical monocytes"   ,  "Intermediate monocytes"   ,   "Classical monocytes",        
                       "mDC",                "pDC" , "CD56+ classical monocytes" )

colnames(signmat) <- colnames(corrmat)
colnames(fdrmat) <- colnames(corrmat)
rownames(signmat) <- rownames(corrmat)
rownames(fdrmat) <- rownames(corrmat)

## correlation map 
heatcolors1 <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(25)
pdf("~/Ped-Covid Studie/FACS/Plots/correlation_ex_COV006.pdf", onefile = T, paper="a4r",width=30, height=20)
corrplot(
  corr = corrmat,
  method="square",
  p.mat = signmat,
  sig.level = c(.001, .01, .05),
  insig = "label_sig",
  pch.cex = .6,
  pch.col = "white",
  col = heatcolors1,
  tl.col = "black",
  tl.cex = 0.75,
  tl.srt = 60,
  cl.pos='n',
  addgrid.col = ifelse(fdrmat < 0.1, "black", rgb(0.6,0.6,0.6,0.2))
)
grid.newpage()
colorlegend(xlim=c(-5, 0), ylim=c(0,5), heatcolors1, c(seq(-1,1,.25)), align="l", vertical=TRUE, addlabels=TRUE)
dev.off()
recordPlot

#plot 
head(melted_cormat)
tiler <- ggplot(data = melted_cormat, aes(x=Var2.x, y=Var1.x, fill=value)) + 
  geom_tile()+
  geom_text(aes(label=stars), color="black", size=4) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" )+
  labs(x = "Cell populations",
       y= "Clinical measurments",
       fill = "Spearman's rho")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust=1, size = 7))

tiler
heatmather <- pheatmap(corrmat,
         cluster_rows = FALSE)


#internal correlations 

#helper functions 
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

cell_cormat <- result$r[1:59,1:59]
upper_tri_cells <- get_upper_tri(cell_cormat)
melted_cormat_cells <- melt(upper_tri_cells, na.rm = TRUE)
melted_cormat_cells$merge <- paste(melted_cormat_cells$Var1, melted_cormat_cells$Var2)

cell_signmat <- result$P[1:59,1:59]
upper_tri_sign <- get_upper_tri(cell_signmat)
melted_cormat_sign <- melt(upper_tri_sign, na.rm = TRUE)
melted_cormat_sign$merge <- paste(melted_cormat_sign$Var1, melted_cormat_sign$Var2)
colnames(melted_cormat_sign)[3] <- "P"

melted_cormat_cells <- merge(melted_cormat_cells, melted_cormat_sign, by = "merge")
melted_cormat_cells <- melted_cormat_cells[,c(1:4,7)]
melted_cormat_cells$stars <- cut(melted_cormat_cells$P, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

kendrick_cellmar <- ggplot(data = melted_cormat_cells, aes(x=Var1.x, y=Var2.x, fill=value)) + 
  geom_tile()+
  geom_text(aes(label=stars), color="black", size=4) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
                       high = "red", space = "Lab" )+
  labs(x = NULL,
       y= NULL,
       fill = "Spearman's rho")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust=1, size = 7),
        legend.position = NULL) +
  scale_x_discrete(position = "top") 


dim(result$r)  

lab_cormat <- result$r[60:86,60:86]
lab_cormat <- lab_cormat[c(1:8,10:15,20,23:27),c(1:8,10:15,20,23:27)]

lower_tri_lab <- get_lower_tri(lab_cormat)
melted_cormat_lab <- melt(lower_tri_lab, na.rm = TRUE)
melted_cormat_lab$merge <- paste(melted_cormat_lab$Var1, melted_cormat_lab$Var2)

lab_signmat <- result$P[60:86,60:86]
lab_signmat <- lab_signmat[c(1:8,10:15,20,23:27),c(1:8,10:15,20,23:27)]

lower_tri_sign <- get_lower_tri(lab_signmat)
melted_cormat_sign_lab <- melt(lower_tri_sign, na.rm = TRUE)
melted_cormat_sign_lab$merge <- paste(melted_cormat_sign_lab$Var1, melted_cormat_sign_lab$Var2)
colnames(melted_cormat_sign_lab)[3] <- "P"

melted_cormat_lab <- merge(melted_cormat_lab, melted_cormat_sign_lab, by = "merge")
melted_cormat_lab <- melted_cormat_lab[,c(1:4,7)]
melted_cormat_lab$stars <- cut(melted_cormat_lab$P, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

delalab <- ggplot(data = melted_cormat_lab, aes(x=Var1.x, y=Var2.x, fill=value)) + 
  geom_tile()+
  geom_text(aes(label=stars), color="black", size=4) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
                       high = "red", space = "Lab" )+
  labs(x = NULL,
       y= NULL,
       fill = "Spearman's rho")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust=1, size = 8),
        legend.position = NULL) +
  scale_y_discrete(position = "right")

delalab



pdf("~/Päd-Covid Studie/FACS/Plots/correlation_matrices_cleaned_new.pdf", onefile = T, paper="a4r",width=20, height=10)

corrplot(
  corr = corrmat,
  method="square",
  p.mat = signmat,
  sig.level = c(.001, .01, .05),
  insig = "label_sig",
  pch.cex = .6,
  pch.col = "white",
  col = heatcolors1,
  tl.col = "black",
  tl.cex = 0.75,
  tl.srt = 60,
  cl.pos='n',
  addgrid.col = ifelse(fdrmat < 0.1, "black", rgb(0.6,0.6,0.6,0.2))
)
grid.newpage()
print(tiler)
grid.newpage()
print(heatmather)
print(kendrick_cellmar)
print(delalab)
dev.off()


lm_eqn <- function(df, y, formula, x){
  m <- lm(as.formula(formula), df);
  eq <- substitute(italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


combined$IL6 <- combined$`IL-6 (<5,9 pg/ml)`
combined$mono_mfi_hladr <- combined$`monocytes_MFI_HLA-DR`
combined$crp <- combined$`CRP mg/dl`
combined$fibr <- combined$`Fibrinogen mg/dl`
combined$krea <- combined$`Kreatinin mg/dl`
combined$antith <- combined$`Antithrombin %`

ggplot(combined, aes_string(x="IL6", y="CM_CD4")) + 
  geom_point(shape = 21,size = 5, fill = "#66C2A5", stroke = 1.6) +
  geom_smooth(color="black", method = "glm", formula = y~x,se = F) + 
  geom_text(aes(label = rownames(combined)))+
  theme_classic() +
  ylim(0,65)+
  geom_text(x = 150, y = 40,label = lm_eqn(combined , CM_CD4 ,"CM_CD4 ~ IL6", IL6), parse = TRUE, size = 8)+
  #labs(y = "Intermediate Monocytes (%)",
   #    x = "Fibrinogen (mg/dl)") +
  theme(   axis.text = element_text( size = 20 ),
           axis.title = element_text( size = 21),
           axis.title.x = element_text(vjust = -0.5),
           axis.title.y = element_text(vjust = 1.2),
           legend.title = element_text(size = 14),
           legend.text = element_text(size = 12)) 


ggsave("intermed_fibr.pdf", height = 6, width = 7)

ggplot(combined, aes_string(x="antith", y="class_mono")) + 
  geom_point(shape = 21,size = 5, fill = "#66C2A5", stroke = 1.6) +
  geom_smooth(color="black", method = "glm", formula = y~x,se = F) + 
  theme_classic() +
  #ylim(0,65)+
  geom_text(x = 75, y = 80,label = lm_eqn(combined , class_mono,"class_mono ~ antith", antith), parse = TRUE, size = 8)+
  #geom_text(aes(label = rownames(combined)))
  labs(y = "Classical Monocytes (%)",
       x = "Antithrombin (%)") +
  theme(   axis.text = element_text( size = 20 ),
           axis.title = element_text( size = 21),
           axis.title.x = element_text(vjust = -0.5),
           axis.title.y = element_text(vjust = 1.2),
           legend.title = element_text(size = 14),
           legend.text = element_text(size = 12)) 

ggsave("class_mono_antithrombin.pdf", height = 6, width = 7)


ex_Cov6 <- combined[-2,]
ggplot(ex_Cov6, aes_string(x="IL6", y="activ_CD8m")) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  ylim(0,65)+
  geom_text(aes(label=rownames(ex_Cov6)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(   axis.text = element_text( size = 14),
           axis.title = element_text( size = 16),
           axis.title.x = element_text(vjust = -0.5),
           axis.title.y = element_text(vjust = 1.2),
           legend.title = element_text(size = 14),
           legend.text = element_text(size = 12)) 

combined$`Fibrinogen mg/dl`
combiend$mono
ggplot(combined, aes(x=`CRP mg/dl`, y=monocyte_mfi_HLA)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(combined)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))


ggplot(ex_Cov6, aes(x=`CRP mg/dl`, y=)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(ex_Cov6)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))


ggplot(combined, aes(x=`CRP mg/dl`, y=CD4)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(combined)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))


ggplot(ex_Cov6, aes(x=`CRP mg/dl`, y=CD4)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(ex_Cov6)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))


ggplot(combined, aes(x=`Fibrinogen mg/dl`, y=intermed_mono)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(combined)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))


ggplot(ex_Cov6, aes(x=`Fibrinogen mg/dl`, y=intermed_mono)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(ex_Cov6)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))


colnames(combined) <- str_replace(colnames(combined), "-", "_")
combined$`Antithrombin %`

ggplot(combined, aes(x=`Antithrombin %`, y=activ_CD4m_mfi_HLA_DR)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(combined)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))


ggplot(ex_Cov6, aes(x=`Antithrombin %`, y=activ_CD4m_mfi_HLA_DR)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(ex_Cov6)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))


ggplot(combined, aes(x=`Antithrombin %`, y=Th2)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(combined)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))


ggplot(ex_Cov6, aes(x=`Antithrombin %`, y=Th2)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(ex_Cov6)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))

ggplot(combined, aes(x=`Thrombozyten G/l`, y=Th2)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(combined)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))


ggplot(ex_Cov6, aes(x=`Thrombozyten G/l`, y=Th2)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(ex_Cov6)),hjust=0, vjust=0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))

combined$pct <- combined$`Procalcitonin ng/ml`
combined$ptt <- combined$`pTT sec`



lm_eqn <- function(df, y, formula, x){
  m <- lm(as.formula(formula), df);
  eq <- substitute(italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}



ggplot(combined, aes(y=class_switch_CD19m, x=ptt)) + 
  geom_point(size = 5) +
  geom_smooth(method = "glm", formula = y~x,se = FALSE) + 
  theme_classic() +
  geom_text(aes(label=rownames(combined)),hjust=0, vjust=0) +
  geom_text(x = 35, y = 15,label = lm_eqn(combined , class_switch_CD19m,"class_switch_CD19m ~ ptt", ptt), parse = TRUE, size = 5)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16, face = "bold" ),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12))





# correlate using pairwise complete observations 
result<-rcorr(as.matrix(ex_Cov6))
result$n
result$r[result$n<6]<-NA # ignore less than 7 observations
rownames(result$r)
result$P %>% view()

#extract relevant columns & rows 
corrmat <- result$r[lab_columns,1:59]

#add statistical significance levels 
melted_cormat <- melt(corrmat)
melted_cormat$merge <- paste(melted_cormat$Var1, melted_cormat$Var2, sep = "_")

signmat <- result$P[lab_columns,1:59]
melted_signmat <-  melt(signmat)
melted_signmat$merge <- paste(melted_signmat$Var1, melted_signmat$Var2, sep = "_")
colnames(melted_signmat)[3] <- "P"

melted_cormat <- merge(melted_cormat, melted_signmat, by = "merge")
melted_cormat <- melted_cormat[,c(1:4,7)]
melted_cormat$stars <- cut(melted_cormat$P, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

#plot 
head(melted_cormat)
tiler <- ggplot(data = melted_cormat, aes(x=Var2.x, y=Var1.x, fill=value)) + 
  geom_tile()+
  geom_text(aes(label=stars), color="black", size=4) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
                       high = "red", space = "Lab" )+
  labs(x = "Cell populations",
       y= "Clinical measurments",
       fill = "Pearson's r")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust=1, size = 7))

tiler
heatmather <- pheatmap(corrmat,
                       cluster_rows = FALSE)

pdf("~/Päd-Covid Studie/FACS/Plots/correlation_matrices_cleaned_EXCOV6.pdf", onefile = T, paper="a4r",width=20, height=10)
print(tiler)
grid.newpage()
print(heatmather)
dev.off()

?rcorr


# correlate using pairwise complete observations 
result<-rcorr(as.matrix(combined), type="spearman")
result$n
result$r[result$n<7]<-NA # ignore less than 7 observations
rownames(result$r)
result$P %>% view()

#extract relevant columns & rows 
corrmat <- result$r[lab_columns,1:59]

#add statistical significance levels 
melted_cormat <- melt(corrmat)
melted_cormat$merge <- paste(melted_cormat$Var1, melted_cormat$Var2, sep = "_")

signmat <- result$P[lab_columns,1:59]
melted_signmat <-  melt(signmat)
melted_signmat$merge <- paste(melted_signmat$Var1, melted_signmat$Var2, sep = "_")
colnames(melted_signmat)[3] <- "P"

melted_cormat <- merge(melted_cormat, melted_signmat, by = "merge")
melted_cormat <- melted_cormat[,c(1:4,7)]
melted_cormat$stars <- cut(melted_cormat$P, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

#plot 
head(melted_cormat)
tiler <- ggplot(data = melted_cormat, aes(x=Var2.x, y=Var1.x, fill=value)) + 
  geom_tile()+
  geom_text(aes(label=stars), color="black", size=4) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
                       high = "red", space = "Lab" )+
  labs(x = "Cell populations",
       y= "Clinical measurments",
       fill = "Spearman")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust=1, size = 7))

tiler
heatmather <- pheatmap(corrmat,
                       cluster_rows = FALSE)

pdf("~/Päd-Covid Studie/FACS/Plots/correlation_matrices_cleaned_spearman.pdf", onefile = T, paper="a4r",width=20, height=10)
print(tiler)
grid.newpage()
print(heatmather)
dev.off()
