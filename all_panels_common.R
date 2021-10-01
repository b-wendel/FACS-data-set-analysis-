library(tidyverse)
library(ggplot2)
library(ggpubr)
library(trelliscopejs)
library(plotly)
library(ggforce)
library(readr)
library(rebus)
library(stringr)
library(hrbrthemes)
library(ggtext)
#library(PCAtools)
library(rstatix)
library(gridExtra)
library(viridis)
library(grid)
#install.packages("pheatmap")
#library(pheatmap)
#install.packages("GGally")
library(GGally)
library(RColorBrewer)
#library(ggbeeswarm)
library(reshape2)
library(fastmap)

setwd("~/Ped-Covid Studie")
load("~/Ped-Covid Studie/FACS/scripts/all_panels_common_WS.RData")
load("~/Ped-Covid Studie/FACS/scripts/20210602_MonoDCBcells.RData")
 
common <- panel1_cov$analysis.identifier.x[panel1_cov$analysis.identifier.x %in% panel2_cov$analysis.identifier.x & 
                                             panel1_cov$analysis.identifier.x %in% P3_cov_complete$analysis.identifier.x  ]
# COV46 has not B-cells panel and gets removed
pall <- panel1_cov[panel1_cov$analysis.identifier.x %in% common,] %>% 
  select(colnames(panel1_cov)[c(6:10,12:13,15,17:19,21,22,25:41,52,53,55,57,58,60,61,63,64,66,67,76,77,79,80,172)]) %>% 
  left_join(panel2_cov[panel2_cov$analysis.identifier.x %in% common,c(6:14,98)], by = "merge_ID") %>%
  left_join(P3_cov_complete[P3_cov_complete$analysis.identifier.x %in% common,c(4:19,24:29,33,116)], by = "merge_ID") %>%
  left_join(mono_dc_bcells[mono_dc_bcells$analysis.identifier.x %in% common,c(2:16, 99)], by = "merge_ID") %>% 
  #metadata
  left_join(P3_cov_complete[P3_cov_complete$analysis.identifier.x %in% common,c(34:125)], by = "merge_ID") %>%
  #add columns
  mutate(Disease_course = cut(max_WHO_classification, breaks=c(-Inf, 0.9, 1.9, 3.9, 5.9,Inf), 
                              label=c("Healthy", "Asymptomatic", "Mild", "Moderate","Severe"))) %>% 
  mutate(Disease_course = factor(Disease_course, levels = c("Healthy","Asymptomatic","Mild","Moderate","Severe"))) %>% 
  mutate(Age_group = cut(Age, breaks=c(-Inf, 0.9, 10.9, 18.9,Inf), 
                         label=c("<1", "1-10", "11-18", ">18"))) %>% 
  mutate(Age_group = factor(Age_group, levels = c("<1","1-10","11-18",">18"))) %>% 
  mutate(Status = factor(Status , levels=c("Healthy", "COVID", "Non-COVID","MISC","Long-COVID"))) %>%
  mutate(sex = ifelse(pall$`Sex_(f=0;_m=1)`==0,"f","m"))

# map status to shorter labels    
m = fastmap()
m$mset("COVID" = "C", "Non-COVID" = "NC", "Healthy" = "H", "MISC" = "M", "Long-COVID" = "LC")

pall$Status[is.na(pall$Status)] <- "Long-COVID"
pall$outlier_code[pall$Status == "Long-COVID"] <- 5

pall$short_status <- unlist(m$mget(as.character(pall$Status)))
pall$short_status <- factor(pall$short_status, levels = c("H", "C", "NC", "MISC","LC"))
  

pall_clean <- pall %>% filter(Age_group != "<1" & outlier_code != 5 & Status != "MISC" )

collapsed <- pall %>% mutate(Age_group = "All_patients")                                                                                                     
collapsed <- rbind(collapsed, pall) 
collapsed$Age_group <- factor(collapsed$Age_group, levels=c("All Patients", "1-10", "11-18", ">18"))

pall_long <- pall_clean %>% gather(key = "parameter", value = "value", c(1:45,47:93))

pall_clean_T0 <- pall_clean %>% filter(day == 0 )
pall_clean_T0 <- pall_clean_T0[-50,]


lm_eqn <- function(df, y, formula, x){
  m <- lm(as.formula(formula), df);
  eq <- substitute(italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


theme_linedraw_nogrid <- theme(panel.grid = element_blank(), 
                               panel.border = element_blank(), 
                               axis.line = element_line(), 
                               strip.background = element_blank(),
                               strip.text = element_text(color="black", face="bold"),
                               plot.title = element_text(hjust=0.5))
theme_linedraw_nogrid_facet <- theme(panel.grid = element_blank(), 
                                     strip.background = element_blank(),
                                     #strip.background = element_rect(fill = "white"),
                                     strip.text = element_text(color="black", face="bold"),
                                     plot.title = element_text(hjust=0.5))

p.labs <- c("All groups", "1-10 y.o.","11-18 y.o.", ">18 y.o." )
names(p.labs) <- c("All Patients", "1-10", "11-18", ">18")
p.labs.truncated <-  c( "1-10 y.o.","11-18 y.o.", ">18 y.o." )
names(p.labs.truncated) <- c("1-10", "11-18", ">18")

#pca analysis ############
#only day 0, no outliers 
pall_norm_pca <- pall %>% filter(outlier_code < 5 & day == 0 )
#only 1-18
pall_norm_pca<- pall_norm_pca %>% filter(Age_group == "1-10" | Age_group == "11-18")

pall_norm_pca <- as.data.frame(pall_norm_pca[,1:78])

pall_norm_pca$merge_ID

rownames(pall_norm_pca) <- pall_norm_pca$merge_ID
pall_norm_pca <- pall_norm_pca[,-46]
rownames(pall_norm_pca)
pall_norm_pca <- as.data.frame(t(as.matrix(pall_norm_pca)))

## ex COV040 wegen fehlender Messungen im Mono-gate  
colnames(pall_norm_pca)
pall_norm_pca <- pall_norm_pca[, -27]

#normalize into z-scores
for (i in 1:length(rownames(pall_norm_pca))) {
  m <- sum(pall_norm_pca[i,])/length(pall_norm_pca[i,])
  std <- sd(pall_norm_pca[i,])
  print(m)
  for (j in 1:length(colnames(pall_norm_pca))) {
    pall_norm_pca[i,j] <- (pall_norm_pca[i,j] - m) / std
  }
}

#only cells
pall_norm_pca_cells <- pall_norm_pca[c(1:27,46:70),]
rownames(pall_norm_pca_cells)

#generate meta data and necessary columns 
meta_pca <- as.data.frame(pall[,c(46,79:170)][pall$outlier_code < 5 & pall$day == 0 ,])
meta_pca <- meta_pca %>% filter(Age_group == "1-10"|Age_group == "11-18")
rownames(meta_pca) <- meta_pca$merge_ID
#check if meta is available for all samples in pca 
rownames(meta_pca) %in% colnames(pall_norm_pca_cells)
# exclude COV40 because of missing Monogate measurements
meta_pca <- meta_pca[-27,]

meta_pca$covid <- ifelse(meta_pca$Status == "COVID", 1, 0)
meta_pca$non_covid <- ifelse(meta_pca$Status=="Non-COVID", 1,0)
meta_pca$healthy <- ifelse(meta_pca$Status == "Healthy", 1,0)

#meta_pca$misc <- ifelse(meta_pca$Status == "MISC", 1,0)

meta_pca$BMI <- as.numeric(sub(",", ".", meta_pca$BMI, fixed = TRUE))
meta_pca$days_symptom_start <- meta_pca$sampling_date - meta_pca$date_of_1st_symptoms
meta_pca$days_symptom_start <- as.numeric(str_remove(meta_pca$days_symptom_start, " days"))
#select relevant columns 
meta_pca <- meta_pca[,c(9,10,13,16,18,30,31,33:42,46,61,81,82,87:89,93:97)]
meta_pca$analysis.identifier <- rownames(meta_pca)
#add missing WHO classifications 
meta_pca$WHO_on_date[40:42] <- meta_pca$max_WHO_classification[40:42] 

colnames(pall_norm_pca) == rownames(meta_pca)


#PCA - cells 

pca_norm_cells <- pca(pall_norm_pca_cells, metadata = meta_pca)



scree <- screeplot(pca_norm_cells, 
          title = 'SCREE plot',
          subtitle = 'PC1-12 explain 80% of variance',
          components = getComponents(pca_norm_cells, 1:25),
          hline = 80, vline = 12, axisLabSize = 14, titleLabSize = 20,
          returnPlot = FALSE)

colnames(meta_pca)

eigcor <- eigencorplot(pca_norm_cells, 
             metavars = c(colnames(meta_pca)[c(1:3,6,18,19, 23:24,26:29)]),
             main = 'PC1-10 - clinical correlations')

aggregate(x = meta_pca$Age,                # Specify data column
          by = list(meta_pca$sex),              # Specify group indicator
          FUN = mean)


cor(meta_pca$WHO_on_date[c(1:39,43:47)], meta_pca$`medication(0=none,1=NSAID,2=AB,3=IS,4=multi`[c(1:39,43:47)],method = "pearson", )
cor.test(meta_pca$WHO_on_date[c(1:39,43:47)], meta_pca$`medication(0=none,1=NSAID,2=AB,3=IS,4=multi`[c(1:39,43:47)], method = "pearson")

#biplots PC1 / PC2 
bp_pc12_age <-biplot(pca_norm_cells, 
       x = 'PC1', y = 'PC2',
       colby = 'Age_group',
       #colkey = c('COVID' = '#BC5090', 'Non-COVID' = '#FFA600', 'Healthy' = "#0381A1", "MISC" = "blue"),
       #encircle = TRUE,encircleAlpha = 1/6,
       lab = meta_pca$analysis.identifier,
       legendPosition = 'top', legendLabSize = 12, legendIconSize = 8.0,legendTitleSize = 14,
       shape = 'Status',
       #shapekey = c("mild"=19, 'moderate'=15,"severe" = 17),
       #lab = NULL,
       drawConnectors = FALSE,
       axisLabSize = 16,
       title = 'PC1 / PC2',
       titleLabSize = 18,
       pointSize = 5)

bp_pc12_status <- biplot(pca_norm_cells, 
       x = 'PC1', y = 'PC2',
       colby = 'Status',
       colkey = c('COVID' = '#BC5090', 'Non-COVID' = '#FFA600', 'Healthy' = "#0381A1", "MISC" = "blue"),
       #encircle = TRUE,encircleAlpha = 1/6,
       lab = meta_pca$analysis.identifier,
       legendPosition = 'top', legendLabSize = 12, legendIconSize = 8.0,legendTitleSize = 14,
       shape = 'Disease_course',
       shapekey = c("mild"=19, 'moderate'=15,"severe" = 17),
       drawConnectors = FALSE,
       axisLabSize = 16,
       title = 'PC1 / PC2',
       titleLabSize = 18,
       pointSize = 5)


#biplots of PC's
bp_pc34_1  <- biplot(pca_norm_cells, 
                   x = 'PC3', y = 'PC4',
                   colby = 'Age_group',
                   #colkey = c('COVID' = '#BC5090', 'Non-COVID' = '#FFA600', 'Healthy' = "#0381A1", "MISC" = "blue"),
                   #encircle = TRUE,encircleAlpha = 1/6,
                   legendPosition = 'top', legendLabSize = 12, legendIconSize = 8.0,legendTitleSize = 14,
                   shape = 'Disease_course', shapekey = c("mild"=19, 'moderate'=15,"severe" = 17),
                   lab = meta_pca$analysis.identifier,
                   drawConnectors = FALSE,
                   axisLabSize = 16,
                   title = 'PC3 / PC4',
                   titleLabSize = 18,
                   pointSize = 5)
#subtitle = 'PC1 versus PC2')

bp_pc34_1

bp_pc34_2 <- biplot(pca_norm_cells, 
                  x = 'PC3', y = 'PC4',
                  colby = 'Status',
                  colkey = c('COVID' = '#BC5090', 'Non-COVID' = '#FFA600', 'Healthy' = "#0381A1", "MISC" = "blue"),
                  #encircle = TRUE,encircleAlpha = 1/6,
                  legendPosition = 'top', legendLabSize = 12, legendIconSize = 8.0,legendTitleSize = 14,
                  shape = 'Disease_course', shapekey = c("mild"=19, 'moderate'=15,"severe" = 17),
                  lab = meta_pca$analysis.identifier,
                  drawConnectors = FALSE,
                  axisLabSize = 16,
                  title = 'PC3 / PC4',
                  titleLabSize = 18,
                  pointSize = 5)
#subtitle = 'PC3 versus PC4')

bp_pc34_2




#visualize age dependence
pc1.age <- as.data.frame(cbind(meta_pca$Age,pca_norm_cells$rotated$PC1))
colnames(pc1.age) <- c("Age", "PC1")

lm_eqn <- function(df, pc, variable){
  m <- lm(PC1 ~ Age, df);
  eq <- substitute(italic(pc) == a + b %.% italic(variable)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

 
pc1_age <- ggplot(pc1.age, aes(x=Age, y=PC1)) +
  geom_point() +
  geom_smooth(method = lm) +
  geom_text(x = 10, y = (-6),label = lm_eqn(pc1.age, PC1, Age), parse = TRUE)+
  theme_classic() +
  labs(title = "Age ~ PC1", 
       x = "Age",
       y = "PC1 ( 20.31% )")+
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 14, face = "bold" ))

pc1_age
#top loadings
pc1.loadings <- as.data.frame(cbind(pca_norm_cells$loadings$PC1, pca_norm_cells$xvars))
pc1.loadings <- pc1.loadings[order(pc1.loadings, decreasing=TRUE),][53:104,]
colnames(pc1.loadings) <- c("Loadings", "Celltype")
rownames(pc1.loadings)
PC1_toploadings <- pc1.loadings[c(1:10,27:36),]
PC1_toploadings$Celltype <-factor(PC1_toploadings$Celltype, levels = PC1_toploadings$Celltype)
PC1_toploadings$Loadings <- as.numeric(PC1_toploadings$Loadings)

p_pc1_pl <- ggplot(PC1_toploadings, aes(x=Celltype, y=Loadings)) +
  geom_segment( aes(x=Celltype, xend=Celltype, y=0, yend=Loadings), color="grey")+
  geom_point( color="orange", size=4) +
  geom_hline(yintercept = 0, alpha=0.6) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, face = "bold"),
    axis.text = element_text( size = 11 ),
    axis.title = element_text( size = 16),
    strip.text = element_text( size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  coord_flip() +
  labs(title ="PC1 Loadings", 
       caption ="",
       x = "",
       y = "Loading")


p_pc1_pl

#visualize WHO correlation & Status / PC2
pc2.who <- as.data.frame(cbind(meta_pca$WHO_on_date, pca_norm_cells$rotated$PC2))
pc2.who <- cbind(pc2.who, meta_pca$Status)
colnames(pc2.who) <- c("WHO_classification", "PC2", "Status")
pc2.who$WHO_classification <- as.factor(pc2.who$WHO_classification)

comparisons_pc2 <- list(c("0", "1"),c("0", "2"),c("0", "3"),c("0", "4"),c("0", "5"),c("0", "6"),
                        c("1", "2"),c("1", "3"),c("1", "4"),c("1", "5"),c("1", "6"),c("2", "3"),
                        c("2", "4"),c("2", "5"),c("2", "6"),c("3", "4"),c("3", "5"),c("3", "6"),
                        c("4", "5"),c("4", "6"),c("5", "6"))
 
pc2_who <- ggplot(pc2.who %>% filter(is.na(WHO_classification) == FALSE), aes(x=WHO_classification, y=PC2, fill = WHO_classification)) +
  geom_boxplot() +
  stat_compare_means(comparisons = comparisons_pc2, aes(label =..p.format..),
                     paired = FALSE) +
  theme_classic() +
  labs(title = "WHO classification ~ PC2", 
       x = "WHO grade",
       y = "PC2 ( 11.34% )")+
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.title = element_text( size =14, face = "bold" ))

pc2_who

my_comparisons <- list( c("COVID", "Non-COVID"), c("COVID", "Healthy"), c("Non-COVID", "Healthy"))


pc2_status <- ggplot(pc2.who, aes(x=Status, y=PC2, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, aes(label =..p.format..),
                     method = "wilcox.test" , paired = FALSE) +
  theme_classic() +
  scale_fill_manual(labels = c("Healthy", "Non-Covid","COVID"),
  values = c("#0381A1", "#FFA600","#BC5090"))+
  labs(title = "Status ~ PC2", 
       x = "Status",
       y = "PC2 ( 11.68% )")+
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 10 ),
         axis.title = element_text( size =14, face = "bold" ))

pc2_status

#top loadings pc2
pc2.loadings <- as.data.frame(cbind(pca_norm_cells$loadings$PC2, pca_norm_cells$xvars))
pc2.loadings <- pc2.loadings[order(pc2.loadings, decreasing=TRUE),][53:104,]
colnames(pc2.loadings) <- c("Loadings", "Celltype")
rownames(pc2.loadings)
PC2_toploadings <- (pc2.loadings[c(1:10,25:34),])

PC2_toploadings$Celltype <-factor(PC2_toploadings$Celltype, levels = PC2_toploadings$Celltype)
PC2_toploadings$Loadings <- as.numeric(PC2_toploadings$Loadings)

pc2_pl <- ggplot(PC2_toploadings, aes(x=Celltype, y=Loadings)) +
  geom_segment( aes(x=Celltype, xend=Celltype, y=0, yend=Loadings), color="grey")+
  geom_point( color="orange", size=4) +
  geom_hline(yintercept = 0, alpha=0.6) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, face = "bold"),
    axis.text = element_text( size = 11 ),
    axis.title = element_text( size = 16),
    strip.text = element_text( size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  coord_flip() +
  labs(title ="PC2 Loadings", 
       caption ="",
       x = "",
       y = "Loading")

pc2_pl

#COVID correlation 
pc3.covid <- as.data.frame(cbind(as.factor(meta_pca$WHO_on_date), pca_norm_cells$rotated$PC3))
pc3.covid <- cbind(pc3.covid, meta_pca$Status)
colnames(pc3.covid) <- c("WHO", "PC3","Status")
pc3.covid$Status<- as.factor(pc3.covid$Status)

pc3_status <- ggplot(pc3.covid, aes(x=Status, y=PC3, fill =Status)) + 
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test" , comparisons = my_comparisons, aes(label =..p.format..)) +
  theme_classic() +
  scale_fill_manual(labels = c("Healthy", "Non-Covid","COVID","MISC"),
                    values = c("#0381A1", "#FFA600","#BC5090","blue"))+
  labs(title = "Status ~ PC3", 
       x = "Status",
       y = "PC3 ( 7.54% )")+
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 10 ),
         axis.title = element_text( size =14, face = "bold" ))

pc3_status

#top loadings
pc3.loadings <- as.data.frame(cbind(pca_norm_cells$loadings$PC3, pca_norm_cells$xvars))
pc3.loadings <- pc3.loadings[order(pc3.loadings$V1, decreasing=TRUE),]
colnames(pc3.loadings) <- c("Loadings", "Celltype")
rownames(pc3.loadings)
pc3_toploadings <- pc3.loadings[c(1:10,25:34),]

pc3_toploadings$Celltype <-factor(pc3_toploadings$Celltype, levels = pc3_toploadings$Celltype)
pc3_toploadings$Loadings <- as.numeric(pc3_toploadings$Loadings)

pc3_pl <- ggplot(pc3_toploadings, aes(x=Celltype, y=Loadings)) +
  geom_segment( aes(x=Celltype, xend=Celltype, y=0, yend=Loadings), color="grey")+
  geom_point( color="orange", size=4) +
  geom_hline(yintercept = 0, alpha=0.6) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, face = "bold"),
    axis.text = element_text( size = 11 ),
    axis.title = element_text( size = 16),
    strip.text = element_text( size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  coord_flip() +
  labs(title ="PC3 Loadings", 
       caption ="",
       x = "",
       y = "Loading")

pc3_pl
# pc3

pc4.covid <- as.data.frame(cbind(as.factor(meta_pca$WHO_on_date), pca_norm_cells$rotated$PC4))
pc4.covid <- cbind(pc4.covid, meta_pca$Status)
colnames(pc4.covid) <- c("WHO", "PC4","Status")
pc4.covid$Status<- as.factor(pc4.covid$Status)

pc4_status <- ggplot(pc4.covid, aes(x=Status, y=PC4, fill =Status)) + 
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test" , comparisons = my_comparisons, aes(label =..p.format..)) +
  theme_classic() +
  scale_fill_manual(labels = c("Healthy", "Non-Covid","COVID","MISC"),
                    values = c("#0381A1", "#FFA600","#BC5090","blue"))+
  labs(title = "Status ~ PC4", 
       x = "Status",
       y = "PC4 ( 7.09% )")+
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 10 ),
         axis.title = element_text( size =14, face = "bold" ))

pc4_status
#top loadings
pc4.loadings <- as.data.frame(cbind(pca_norm_cells$loadings$PC4, pca_norm_cells$xvars))
pc4.loadings <- pc4.loadings[order(pc4.loadings$V1, decreasing=TRUE),]
colnames(pc4.loadings) <- c("Loadings", "Celltype")
rownames(pc4.loadings)
pc4_toploadings <- pc4.loadings[c(1:10,33:42),]

pc4_toploadings$Celltype <-factor(pc4_toploadings$Celltype, levels = pc4_toploadings$Celltype)
pc4_toploadings$Loadings <- as.numeric(pc4_toploadings$Loadings)

pc4_pl <- ggplot(pc4_toploadings, aes(x=Celltype, y=Loadings)) +
  geom_segment( aes(x=Celltype, xend=Celltype, y=0, yend=Loadings), color="grey")+
  geom_point( color="orange", size=4) +
  geom_hline(yintercept = 0, alpha=0.6) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, face = "bold"),
    axis.text = element_text( size = 11 ),
    axis.title = element_text( size = 16),
    strip.text = element_text( size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  coord_flip() +
  labs(title ="PC4 Loadings", 
       caption ="",
       x = "",
       y = "Loading")

pc4_pl



# Heatmap #####
my_sample_col <- meta_pca[,c(20,4)]
unique(pall$merge_ID[pall$Status == "Healthy"])
pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "Healthy"])]
pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "Non-COVID"])]
pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "COVID"])]
pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "MISC"])]

pall_norm_ordered <- pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "Healthy" & pall$Age_group == "1-10"])] %>%
  cbind(pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "Healthy" & pall$Age_group == "11-18"])]) %>% 
  cbind(pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "Non-COVID" & pall$Age_group == "1-10"])])  %>% 
  cbind(pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "Non-COVID" & pall$Age_group == "11-18"])]) %>% 
  cbind(pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "COVID" & pall$Age_group == "1-10"])]) %>% 
  cbind(pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "COVID" & pall$Age_group == "11-18"])])  %>% 
  cbind(pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "MISC" & pall$Age_group == "1-10"])])

pall_norm_ordered_age <- pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "Healthy" & pall$Age_group == "1-10"])] %>%
  cbind(pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "Non-COVID" & pall$Age_group == "1-10"])])  %>% 
  cbind(pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "COVID" & pall$Age_group == "1-10"])]) %>% 
  cbind(pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "MISC" & pall$Age_group == "1-10"])]) %>% 
  cbind(pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "Healthy" & pall$Age_group == "11-18"])]) %>% 
  cbind(pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "Non-COVID" & pall$Age_group == "11-18"])]) %>% 
  cbind(pall_norm_pca_cells[,colnames(pall_norm_pca_cells) %in% unique(pall$merge_ID[pall$Status == "COVID" & pall$Age_group == "11-18"])])  
  
my_colour = list(
  Status = c("COVID" = "#BC5090", "Non-COVID" ="#FFA600", "Healthy" = "#0381A1", "MISC" = "blue"),
  Age_group = c("1-10" = "#CC5500", "11-18" = "lightblue"))


heat_colors <- brewer.pal(6, "YlOrRd")

hm <- pheatmap(pall_norm_ordered,
               annotation_col = my_sample_col,
               cluster_cols=TRUE, 
               cluster_rows = TRUE,
               cutree_cols = 6,
               annotation_colors = my_colour)

hm_age <-pheatmap(pall_norm_ordered_age ,
                 annotation_col = my_sample_col,cluster_cols=FALSE, cluster_rows = TRUE,
                  annotation_colors = my_colour)

#t-test for all populations & waterfall #####

pall.ttest.covid.healthy <- pall_long %>% filter(Status == "Healthy" | Status == "COVID" ) %>%
  filter(day == 0 & outlier_code < 5) %>%
  filter(Age_group == "1-10" | Age_group =="11-18") %>%
  group_by(parameter)%>%
  t_test(value ~ Status_ordered)%>%
  adjust_pvalue(method = "BH") %>%
  add_significance()


#waterfall of all populations + MFI of all age groups 
st_wf <- pall.ttest.covid.healthy[,c(1,7,9:11)] %>% filter(p.adj < .1)
st_wf <- st_wf[order(st_wf$statistic),]
st_wf$parameter <- as.factor(st_wf$parameter)
st_wf$start <- 0

pattern_mfi <- or("MFI", "mfi")
st_wf_exMFI <- st_wf %>% filter(!str_detect(parameter, pattern_mfi))
st_wf_exMFI$id <- c(1:16)
st_wf_exMFI$parameter <- str_replace_all(st_wf_exMFI$parameter, "_"," ")
st_wf_exMFI$parameter <- str_replace_all(st_wf_exMFI$parameter, "pos","+")

st_wf_MFI <- st_wf %>% filter(str_detect(parameter, pattern_mfi))
st_wf_MFI$id <- c(1:10)
st_wf_MFI$parameter <- str_replace_all(st_wf_MFI$parameter, "_"," ")
st_wf_MFI$parameter <- str_replace_all(st_wf_MFI$parameter, "pos","+")
st_wf_MFI$parameter <- str_replace(st_wf_MFI$parameter, "mfi", "/")
st_wf_MFI$parameter <- str_replace(st_wf_MFI$parameter, "MFI", "/")

st_wf$id <- c(1:26)
st_wf$parameter <- str_replace_all(st_wf$parameter, "_"," ")
st_wf$parameter <- str_replace_all(st_wf$parameter, "pos","+")
st_wf$parameter <- str_replace(st_wf$parameter, "mfi", "/")
st_wf$parameter <- str_replace(st_wf$parameter, "MFI", "/")



waterfall_covid_healthy <- st_wf %>% mutate(parameter = forcats::fct_inorder(parameter)) %>%
ggplot() +
  geom_rect(aes(x = parameter, xmin = id -0.45, xmax = id+0.45,ymin = start, ymax= statistic,fill = statistic)) +
  geom_hline(yintercept = 0, alpha=0.8) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
  theme_classic()+
  scale_fill_gradient2()+
  labs(title ="Distinguishing features - COVID vs Healthy", 
       #subtitle = "all panels - all markers - Age: 1-18", 
       caption ="FDR < 0.1",
       x = "Feature",
       y = "t-statistic") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust=1, size = 9),
        plot.caption = element_text(size = 14))

waterfall_covid_healthy_cells <- st_wf_exMFI%>% mutate(parameter = forcats::fct_inorder(parameter)) %>%
  ggplot() +
  geom_rect(aes(x = parameter, xmin = id -0.45, xmax = id+0.45,ymin = start, ymax= statistic,fill = statistic)) +
  geom_hline(yintercept = 0, alpha=0.8) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
  theme_classic()+
  scale_fill_gradient2()+
  labs(title ="Distinguishing features - COVID vs Healthy", 
       #subtitle = "all panels - all markers - Age: 1-18", 
       caption ="FDR < 0.1",
       x = "Feature",
       y = "t-statistic") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust=1, size = 9),
        plot.caption = element_text(size = 14))

waterfall_covid_healthy_MFI <- st_wf_MFI %>% mutate(parameter = forcats::fct_inorder(parameter)) %>%
  ggplot() +
  geom_rect(aes(x = parameter, xmin = id -0.45, xmax = id+0.45,ymin = start, ymax= statistic,fill = statistic)) +
  geom_hline(yintercept = 0, alpha=0.8) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
  theme_classic()+
  scale_fill_gradient2()+
  labs(title ="Distinguishing features - COVID vs Healthy", 
       #subtitle = "all panels - all markers - Age: 1-18", 
       caption ="FDR < 0.1",
       x = "Feature",
       y = "t-statistic") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust=1, size = 9),
        plot.caption = element_text(size = 14))

waterfall_covid_healthy_cells
waterfall_covid_healthy_MFI

png(filename="waterfall_covid_healthy2.png", width=600, height=320)
print(waterfall_covid_healthy)
dev.off()

png(filename="waterfall_covid_healthy_cells.png", width=600, height=320)
print(waterfall_covid_healthy_cells)
dev.off()

png(filename="waterfall_covid_healthy_MFI.png", width=600, height=320)
print(waterfall_covid_healthy_MFI)
dev.off()



# correlation alter ####
pall.correlation.healthy <- pall_long %>% filter(Status == "Healthy") %>%
  filter(Age_group == "1-10" | Age_group =="11-18") %>%
  group_by(parameter)%>%
  cor_test(value, Age, method="spearman") %>% 
  filter(p < .2) %>% 
  filter()

write_xlsx(pall.correlation.healthy, "healthy_spearman_age_populations_unfiltered.xlsx")

ggplot(pall %>%filter(Status == "Healthy" | Status == "COVID") %>%
         filter(Age_group == "1-10" | Age_group =="11-18") %>% filter(day == 0), aes(x= Age, y= naive_CD4, color = Status)) +
  geom_point() +
  geom_text (aes (label = analysis.identifier.x))
 

# stacking #####
ggplot(pall_long %>% filter(day == 0 & outlier_code < 5) %>%
                filter(Age_group != "<1") %>%
                filter(parameter == "naive_CD4" | parameter == "EM_CD4" | 
                              parameter == "CM_CD4" |parameter == "TEM_CD4"), 
       aes(fill=parameter, y=value, x=Status)) + 
  facet_grid(.~Age_group,
             labeller = labeller(Age_group = p.labs.truncated)) +
  geom_bar(position="fill",  color = "black", stat = "summary")  + 
  theme_classic() +
  labs(title = "CD4 T-cells",
       y= "% / CD4+"
  ) +
  scale_fill_manual(labels = c("CM","EM","Naive","TEM"),
                    values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")) +
  theme( plot.title = element_text(size = 26, face = "bold"),
         axis.text = element_text( size = 20 ),
         axis.text.x = element_blank(),
         axis.title = element_text( size = 21),
         axis.title.x = element_blank(),
         axis.ticks.x = element_blank(),
         strip.text = element_text( size = 23, face="bold"),
         strip.background = element_blank(),
         legend.title = element_blank(),
         legend.text = element_text(size = 18))

ggsave("CD4_stacked_incl_naive.pdf", height = 5, width = 7)


ggplot(pall_long %>%  filter(day == 0 & outlier_code < 5) %>%
                filter(Age_group != "<1") %>%
                filter(parameter == "naive_CD8" | parameter == "EM_CD8" | 
                              parameter == "CM_CD8" |parameter == "TEM_CD8"),
        aes(fill=parameter, y=value, x=Status)) + 
  facet_grid(.~Age_group) +
  geom_bar(position="fill",  color = "black", stat = "summary")  + 
  theme_classic() +
  labs(title = "CD8 T-cells",
       y= "% / CD8+"
       ) +
  scale_fill_manual(labels = c("CM","EM","Naive","TEM"),
                    values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")) +
  theme( plot.title = element_text(size = 26, face = "bold"),
         axis.text = element_text( size = 20 ),
         axis.text.x = element_text(angle = 45, vjust = 1.05, hjust=1),
         axis.title = element_text( size = 21),
         axis.title.x = element_blank(),
         axis.ticks.x = element_blank(),
         strip.text = element_text( size = 23, face="bold"),
         strip.background = element_blank(),
         legend.title = element_blank(),
         legend.text = element_text(size = 18))

  

ggsave("CD8_stacked_incl_naive.pdf", height = 6, width = 7)



ggplot(pall_long %>%  filter(day == 0 & outlier_code < 5) %>%
         filter(Age_group != "<1") %>%
         filter(parameter == "CD4" | parameter == "CD8"),
       aes(fill=parameter, y=value, x=Status)) + 
  #facet_grid(.~Age_group) +
  geom_bar(position="fill", stat="summary", color = "black")  + 
  scale_fill_manual(labels = c("CD4+", "CD8+"),
                    values = c("#66C2A5", "#FC8D62")) +
  theme_classic() +
  labs(y= "CD4/CD8 ratio") +
  theme( #plot.title = element_text(size = 26, face = "bold"),
         axis.text = element_text( size = 20 ),
         axis.text.x = element_text(angle = 45, vjust = 1.05, hjust=1),
         axis.title = element_text( size = 21),
         axis.title.x = element_blank(),
         axis.ticks.x = element_blank(),
         legend.title = element_blank(),
         legend.text = element_text(size = 18))

ggsave("CD4CD8_ratio.pdf", height = 6, width = 4)



ggplot(pall_long %>%  filter(day == 0 & outlier_code < 5) %>%
         filter(Age_group != "<1") %>%
         filter(parameter == "class_mono" | parameter == "intermed_mono" | 
                  parameter == "non_class_mono"),
       aes(fill=parameter, y=value, x=Status)) + 
  facet_grid(.~Age_group,
             labeller = labeller(Age_group = p.labs.truncated)) +
  geom_bar(position="fill", stat= "summary", color="black")  + 
  theme_classic() +
  scale_fill_manual(name = "",
                    labels = c("Classical","Intermediate","Non-Classical"),
                    values = c("#66C2A5","#FC8D62","#8DA0CB"))+
  labs( y= "% / Monocytes")+  
  theme( plot.title = element_text(size = 26, face = "bold"),
         axis.text = element_text( size = 20 ),
         axis.text.x = element_text(angle = 45, vjust = 1.05, hjust=1),
         axis.title = element_text( size = 21),
         axis.title.x = element_blank(),
         axis.ticks.x = element_blank(),
         strip.text = element_text( size = 23, face="bold"),
         strip.background = element_blank(),
         legend.position = "none",
         legend.title = element_blank(),
         legend.text = element_text(size = 18))

ggsave("Monocytes_subpopulations_1.pdf", height = 6, width = 6)


pall_long_bstack <- pall_long %>% filter(parameter == "naive_CD19" | parameter == "non_switch_CD19m" | 
                                           parameter == "class_switch_CD19m" | parameter == "CD27_IgD_DN") 
pall_long_bstack$parameter <- factor(pall_long_bstack$parameter, levels= c("naive_CD19","CD27_IgD_DN","non_switch_CD19m",
                                                                           "class_switch_CD19m"))

ggplot(pall_long_bstack %>%  filter(day == 0 & outlier_code < 5) %>%
                   filter(Age_group != "<1"),
                 aes(fill=parameter, y=value, x=Status)) + 
  facet_grid(.~Age_group,
             labeller = labeller(Age_group = p.labs.truncated)) +
  geom_bar(position="fill", stat= "summary", color="black")  + 
  theme_classic() +
  #scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(name = "",
                    labels = c("Naive", "CD27-IgD-","NS Memory","CS Memory"),
                    values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")) +
  labs(y= "% / CD19+")+  
  theme( plot.title = element_text(size = 26, face = "bold"),
        axis.text = element_text( size = 20 ),
        axis.text.x = element_text(angle = 45, vjust = 1.05, hjust=1),
        axis.title = element_text( size = 21),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text( size = 23, face="bold"),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


ggsave("B.cells.stacked.pdf", height = 6, width = 8)



xpdf("~/Ped-Covid Studie/FACS/Plots/TBMONO_Stackedbarplots.pdf", onefile = T, paper="a4r",width=20, height=20)
print(cd4)
print(cd8)
print(cd4cd8)
print(monos)
print(bcells)
dev.off()







# ratioplots over WHO etc ######
lm_eqn <- function(df, y, formula, x){
  m <- lm(as.formula(formula), df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

ratioplot1 <- ggplot(pall_clean %>% filter(Status == "COVID" & analysis.identifier.x != "COV006") , aes(x=WHO_on_date)) + 
  geom_point(aes(y=CD4, color = "red")) + 
  geom_smooth(aes(y=CD4, color = "red"),method = "glm", formula = y~x) + 
  geom_text(x = 3, y = 70,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), CD4,"CD4 ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5)+
  geom_point(aes(y= CD8, color= "blue")) + 
  geom_smooth(aes(y=CD8, color= "blue"),method = "glm", formula = y~x) +
  geom_text(x = 2, y = 20,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), CD8,"CD8 ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5)+
  labs(title = "CD4 & CD8 / WHO_on_Date",
       subtitle = "data = only COVID, Age 1-18, no outliers ",
       y = "%  / CD3+",
       caption = "CD4 = blue; CD8 = red") +
  theme_classic()
ratioplot1

ratioplot2 <- ggplot(pall_clean %>% filter(Status == "COVID") , aes(x=WHO_on_date)) +
  geom_point(aes(y=non_class_mono, color = "red")) + 
  geom_smooth(aes(y=non_class_mono, color = "red"),method = "glm", formula = y~x) + 
  geom_text(x = 3, y = -5,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), non_class_mono,
                                         "non_class_mono ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5) +
  geom_point(aes(y= intermed_mono, color= "blue")) + 
  geom_smooth(aes(y=intermed_mono, color= "blue"),method = "glm", formula = y~x) +
  geom_text(x = 2, y = 30,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), intermed_mono,
                                         "intermed_mono ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5)+
  geom_point(aes(y= class_mono, color= "yellow")) + 
  geom_smooth(aes(y=class_mono, color= "yellow"),method = "glm", formula = y~x) +
  geom_text(x = 4, y = 90,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), class_mono,
                                         "class_mono ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5) +
  labs(title = "Monocyte types / WHO_on_Date",
       subtitle = "data = only COVID, Age 1-18, no outliers ",
       y = "%  / monocytes+",
       caption = "Non-classial = green; intermediate = red; classical = blue"
  ) +
  theme_classic()


ggplot(pall_clean %>% filter(Status != "Non-COVID" & day == 0) , aes(x=Age)) +
  geom_point(aes(y=Th17, color = Status )) + 
  geom_smooth(aes(y=Th17 ),color = "red",method = "glm", formula = y~x) + 
  #geom_point(aes(y= intermed_mono, color= "blue")) + 
  #geom_smooth(aes(y=intermed_mono, color= "blue"),method = "glm", formula = y~x) +
  theme_classic()





# printint #####
pdf("~/Ped-Covid Studie/FACS/Plots/thebigdude.pdf", onefile = T, paper="a4r",width=20, height=10)
print(scree)
print(eigcor)
print(bp_pc12_age)
print(bp_pc12_status )
print(bp_pc34_1)
print(bp_pc34_2)
print(pc1_age)
print(p_pc1_pl)
print(pc2_who)
print(pc2_status)
print(pc2_pl)
print(pc3_status)
print(pc3_pl)
print(pc4_status)
print(pc4_pl)

grid.newpage()
print(hm)
grid.newpage()
print(hm_age)

print(waterfall_covid_healthy)

print(cd4)
print(cd8)

print(ratioplot1)
print(ratioplot2)
dev.off()






ratioplot3 <- ggplot(pall_clean %>% filter(Status == "COVID" & analysis.identifier.x != "COV006") , aes(x=WHO_on_date)) +
  geom_point(aes(y=activ_CD4m, color = "red")) + 
  geom_smooth(aes(y=activ_CD4m, color = "red"),method = "glm", formula = y~x) + 
  geom_text(x = 4, y = 0,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), activ_CD4m,
                                         "activ_CD4m ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5) +
  geom_point(aes(y= activ_CD8m, color= "blue")) + 
  geom_smooth(aes(y=activ_CD8m, color= "blue"),method = "glm", formula = y~x) +
  geom_text(x = 3, y = 6,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), activ_CD8m,
                                         "activ_CD8m ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5)+
  labs(title = "activated T-cells / WHO_on_Date",
       subtitle = "data = only COVID, Age 1-18, no outliers ",
       y = "%  / CD4+/CD8+",
       caption = "active CD8+ = red; active CD4+ = blue"
  ) +
  theme_classic()

ratioplot3
pdf("~/Ped-Covid Studie/FACS/Plots/activeTcells.pdf", onefile = T, paper="a4r",width=20, height=10)
print(ratioplot3)
dev.off()


colnames(pall_clean)[72:78] <- c("non_class_mono_mfi_HLA_DR",                                               
                                 "intermed_mono_mfi_HLA_DR" ,                                               
                                 "class_mono_mfi_HLA_DR",                                                   
                                 "proinflam_mono_mfi_HLA_DR",                                               
                                 "mDC_mfi_HLA_DR",                                                          
                                 "pDC_mfi_HLA_DR" ,                                                         
                                 "monocytes_MFI_HLA_DR" )

colnames(pall_clean)

ratioplot4 <- ggplot(pall_clean %>% filter(Status == "COVID") , aes(x=WHO_on_date)) +
  geom_point(aes(y=monocytes_MFI_HLA_DR, color = "red")) + 
  geom_smooth(aes(y=monocytes_MFI_HLA_DR, color = "red"),method = "glm", formula = y~x) + 
  geom_text(x = 4, y = 0,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), monocytes_MFI_HLA_DR,
                                        "monocytes_MFI_HLA_DR ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5) +
  labs(title = "Monocytes HLA-DR / WHO_on_Date",
       subtitle = "data = only COVID, Age 1-18, no outliers ",
       y = "MFI"
       #caption = "active CD8+ = red; active CD4+ = blue"
  ) +
  theme_classic()

ratioplot4


ratioplot5 <- ggplot(pall_clean , aes(x=Age)) +
  geom_point(pall_clean %>% filter(Status == "Healthy"),aes(y=monocytes_MFI_HLA_DR, color = "red")) + 
  geom_smooth(pall_clean %>% filter(Status == "Healthy"),aes(y=monocytes_MFI_HLA_DR, color = "red"),method = "glm", formula = y~x) + 
  geom_text(x = 4, y = 0,label = lm_eqn(pall_clean %>% filter(Status == "Healthy"), monocytes_MFI_HLA_DR,
                                        "monocytes_MFI_HLA_DR ~ Age", Age), parse = TRUE, size = 5) +
  geom_point(pall_clean %>% filter(Status == "COVID"),aes(y=monocytes_MFI_HLA_DR, color = "red")) + 
  geom_smooth(pall_clean %>% filter(Status == "COVID"),aes(y=monocytes_MFI_HLA_DR, color = "red"),method = "glm", formula = y~x) + 
  geom_text(x = 4, y = 0,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), monocytes_MFI_HLA_DR,
                                        "monocytes_MFI_HLA_DR ~ Age", Age), parse = TRUE, size = 5) +
  labs(title = "Monocytes HLA-DR / Age",
       subtitle = "data = only COVID, Age 1-18, no outliers ",
       y = "MFI"
       #caption = "active CD8+ = red; active CD4+ = blue"
  ) +
  theme_classic()
ratioplot5


ratioplot6 <- ggplot(pall_clean %>% filter(Status == "COVID") , aes(x=WHO_on_date)) +
  geom_point(aes(y=CD3, color = "red")) + 
  geom_smooth(aes(y=CD3, color = "red"),method = "glm", formula = y~x) + 
  geom_text(x = 3, y = 30,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), CD3,
                                        "CD3 ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5) +
  geom_point(aes(y=CD19, color = "red")) + 
  geom_smooth(aes(y=CD19, color = "red"),method = "glm", formula = y~x) + 
  geom_text(x = 3, y = 0,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), CD19,
                                        "CD19 ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5) +
  labs(title = "CD3 & CD19 / WHO_on_Date",
       subtitle = "data = only COVID, Age 1-18, no outliers ",
       y = "%  / CD45+"
  ) +
  theme_classic()
ratioplot4




# fever over intermed monocytes #####

pall_clean$fever <- pall_clean$`Fever_ever_appeared_during_period_of_illness_(y=1/n=0)`
pall_clean$fever <- ifelse(pall_clean$fever == 1, "yes", "no")
pall_clean$mono_mfi_HLA_DR <- pall_clean$`monocytes_MFI_HLA-DR`

ggplot(pall_clean %>% filter( Status != "Healthy", day == 0 & 
                               analysis.identifier.x != "COV088" & Age_group != ">18"), 
       aes(x= as.factor(fever), y =mono_mfi_HLA_DR , fill =as.factor(fever)))+ 
  geom_boxplot(width = 0.35,  alpha = 0.4)+
  geom_point(aes(color = Status),#shape = 21,
    size = 5, #stroke = 1.2, 
    alpha = 1,position = position_jitter(width = .15)) +
  stat_compare_means(aes(label = paste("P = ",..p.format.., sep = "")),method="wilcox.test", label.x = 1.5, 
                     label.y = 13000, size = 8) +
  scale_fill_manual(values=c("#8DA0CB","#E78AC3") )+
  #geom_text(aes(label = analysis.identifier.x)) +
  labs(y = "Monocyte HLA-DR expression (MFI)",
       x = "Fever")+
  theme_linedraw() + 
  theme_linedraw_nogrid_facet  +
  theme( plot.title = element_text(size = 26, face = "bold"),
         axis.text = element_text( size = 20 ),
         axis.title = element_text( size = 21),
         legend.position = "right",
         axis.ticks.x = element_blank())

ggsave("mono_hladr_fever.pdf", height = 6, width = 4)

# mono mfis waterfall ####
#exclude COV40 because of missing monogate
pall_clean_covid <- pall_clean_covid[-19,]
mfis <- aggregate(x = pall_clean_covid[,72:74],                # Specify data column
          by = list(pall_clean_covid$Age_group, pall_clean_covid$Status),              # Specify group indicator
          FUN = mean,na.rm=TRUE, na.action=NULL)

se <- function(z) {
  sd(z) / sqrt(length(z))
}

ses <- aggregate(x = pall_clean_covid[,72:74],                # Specify data column
          by = list(pall_clean_covid$Age_group, pall_clean_covid$Status),
          FUN = function(x) se(x))

mfis_u10 <- melt(1- (mfis[4,3:5] / mfis[1,3:5]))
mfis_o10 <- melt(1- (mfis[5,3:5] / mfis[2,3:5]))
mfis_o18 <- melt(1- (mfis[6,3:5] / mfis[3,3:5]))

wf_mfi <- rbind(mfis_u10, mfis_o10,mfis_o18)
wf_mfi$Age_group <- c("COVID 1-10","COVID 1-10","COVID 1-10","COVID 11-18","COVID 11-18","COVID 11-18","COVID >18","COVID >18","COVID >18")
wf_mfi$Age_group <- factor(wf_mfi$Age_group, levels= c("COVID 1-10", "COVID 11-18", "COVID >18"))
wf_mfi$variable <- as.factor(wf_mfi$variable)
wf_mfi$variable <- str_remove(wf_mfi$variable, "_mfi_HLA-DR")
wf_mfi$variable <- str_replace_all(wf_mfi$variable, "_", " ")

wf_mfi$value <- (wf_mfi$value - (2*wf_mfi$value)) *100

wf_mfi$start <- 0
wf_mfi$id <- c(1,2,3,1,2,3,1,2,3)
wf_mfi$direction <- ifelse(wf_mfi$value > 0, "positive", "negative")

wf_mfi %>% mutate(variable = forcats::fct_inorder(variable)) %>%
  ggplot() +
  facet_grid(.~ Age_group)+
  geom_rect(aes(x = variable, xmin = id -0.45, xmax = id+0.45,ymin = start, ymax= value,fill = direction
                ), color = "#252525", alpha = .7 , binwidth = 0.5) +
  geom_hline(yintercept = 0, alpha=0.8) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
  theme_classic() +
  scale_fill_manual(values = c(scales::muted("blue"),
                       scales::muted("red"))) + 
  labs(y = "% change of HLA-DR (MFI)
relative to healthy controls ") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust=1, size = 9),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text  = element_text(size=12, face="bold"),
        plot.caption = element_text(size = 10),
        axis.title.x = element_blank())

waterfall_covid_healthy


## tiles of mono MFIs HLA-DR ####

mono_mfi_z <- as.data.frame(pall_clean_T0[,c("merge_ID","non_class_mono_mfi_HLA-DR",
                                          "intermed_mono_mfi_HLA-DR","class_mono_mfi_HLA-DR")])


rownames(mono_mfi_z) <- mono_mfi_z$merge_ID
# ex COV040 
rownames(mono_mfi_z)
mono_mfi_z <- mono_mfi_z[-31,-1]

#normalize into z-scores
for (i in 1:ncol(mono_mfi_z)) {
  m <- sum(mono_mfi_z[,i])/length(mono_mfi_z[,i])
  std <- sd(mono_mfi_z[,i])
  print(m)
  for (j in 1:nrow(mono_mfi_z)) {
    mono_mfi_z[j,i] <- (mono_mfi_z[j,i] - m) / std
  }
}

colnames(pall_clean)

mono_mfi_z$merge_ID <- rownames(mono_mfi_z)
mono_mfi_z <- mono_mfi_z %>% left_join(pall_clean[,c("merge_ID","Status", "Age_group")],by = "merge_ID")
colnames(mono_mfi_z) <- str_replace(colnames(mono_mfi_z), "-", "_")
colnames(mono_mfi_z)

mean_z_scores <- aggregate(cbind(class_mono_mfi_HLA_DR,intermed_mono_mfi_HLA_DR,
                non_class_mono_mfi_HLA_DR) ~ Age_group + Status, 
          data = mono_mfi_z, FUN = mean, na.rm = TRUE)
mean_z_scores <- mean_z_scores %>% gather(key = "parameter", value = "value", c(3:5))
mean_z_scores$name <- paste(mean_z_scores$Status, mean_z_scores$Age_group)
mean_z_scores$name2 <- paste(mean_z_scores$parameter, mean_z_scores$Age_group)
mean_z_scores$name2 <-factor(mean_z_scores$name2, levels=c("class_mono_mfi_HLA_DR 1-10","intermed_mono_mfi_HLA_DR 1-10",
                                                           "non_class_mono_mfi_HLA_DR 1-10",
                                                           "class_mono_mfi_HLA_DR 11-18","intermed_mono_mfi_HLA_DR 11-18",
                                                           "non_class_mono_mfi_HLA_DR 11-18",
                                                           "class_mono_mfi_HLA_DR >18" , "intermed_mono_mfi_HLA_DR >18",
                                                           "non_class_mono_mfi_HLA_DR >18"))
mean_z_scores$Status <- factor(mean_z_scores$Status, levels=c("Non-COVID", "COVID", "Healthy"))
mean_z_scores$parameter <- str_remove(mean_z_scores$parameter, "_mfi_HLA_DR")
mean_z_scores$parameter <- str_replace(mean_z_scores$parameter, "_", " ")
mean_z_scores$merge <-paste(mean_z_scores$Age_group, mean_z_scores$parameter, mean_z_scores$Status)

pall_long$parameter <- str_replace(pall_long$parameter, "-", "_")
mfis_sign <- pall_long %>% filter(Status == "Healthy" | Status == "COVID" | Status =="Non-COVID") %>%
  filter(day == 0 & outlier_code < 5) %>%
  filter( Age_group != "<1") %>%
  filter(parameter == "class_mono_mfi_HLA_DR"| parameter== "intermed_mono_mfi_HLA_DR"|
                       parameter =="non_class_mono_mfi_HLA_DR") %>% 
  group_by(parameter, Age_group)%>%
  wilcox_test(value  ~ Status)%>%
  adjust_pvalue(method = "BH") %>% 
  filter((group1 == "Healthy" & group2=="Non-COVID")|(group1 == "Healthy" & group2=="COVID")) %>%
  mutate(parameter = str_remove(parameter, "_mfi_HLA_DR")) %>% 
  mutate(parameter = str_replace(parameter, "_", " ")) %>% 
  mutate(merge = paste(Age_group, parameter, group2)) 
mfis_sign <- mfis_sign %>% 
  mutate(stars = cut(mfis_sign$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))) ## cut??

mean_z_scores <- mean_z_scores %>% left_join(mfis_sign[,9:13], by = "merge") 
mean_z_scores$stars[is.na(mean_z_scores$stars)] <- "" 
mean_z_scores$p.adj[is.na(mean_z_scores$p.adj)] <- 1
mean_z_scores$FDR <- as.factor(ifelse(mean_z_scores$p.adj < .1, 1,0))

class(mean_z_scores$FDR)

mean_z_scores$parameter[mean_z_scores$parameter == "class mono"] <- "Classical"
mean_z_scores$parameter[mean_z_scores$parameter == "intermed mono"] <- "Intermediate"
mean_z_scores$parameter[mean_z_scores$parameter == "non class_mono"] <- "Non-Classical"



ggplot(data = mean_z_scores, aes(x=parameter, y=Status, 
                                 fill=value)) + 
  facet_grid(.~Age_group,
             labeller = labeller(Age_group = p.labs.truncated)) +
  geom_tile(color = "white", size = 1.2)+
  geom_tile(aes(color = FDR),fill = NA, size=1.7)+
  coord_equal() +
  geom_text(aes(label=stars), color="white", size=10) +
  scale_fill_gradient2(name = "HLA-DR
expression", #labels = c("1" = "high", "-1" = "low"),
                       breaks = c(1,-0.9), labels = c("high", "low"),
                       midpoint = 0, low = scales::muted("blue"), mid = "gray96",
                       high = scales::muted("red"), space = "Lab", n.breaks = 2)+
  scale_color_manual(name ="FDR", values=c("white","black"),
                     labels = c("0"="n.s","1"="FDR<0.1")) +
 #labs(x = "Cell populations",
   #    y= "Clinical measurments",
    #   fill = "Spearman's rho")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05,
                              hjust=1),
    axis.text = element_text(size = 20),
    strip.background = element_blank(),
    strip.text = element_text(size =23, face="bold"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title =element_blank(),
    legend.position = "right", 
    legend.text = element_text(size=16),
    legend.title = element_text(size=17))

ggsave("Mono_Hla_Dr_tiles_large.pdf", height = 6, width = 10)


# DCs waterfall MFI####
mfis <- aggregate(x = pall_clean_covid[,76:77],                # Specify data column
                  by = list(pall_clean_covid$Age_group, pall_clean_covid$Status),              # Specify group indicator
                  FUN = mean,na.rm=TRUE, na.action=NULL)

se <- function(z) {
  sd(z) / sqrt(length(z))
}

ses <- aggregate(x = pall_clean_covid[,76:77],                # Specify data column
                 by = list(pall_clean_covid$Age_group, pall_clean_covid$Status),
                 FUN = function(x) se(x))

mfis_u10 <- melt(1- (mfis[4,3:4] / mfis[1,3:4]))
mfis_o10 <- melt(1- (mfis[5,3:4] / mfis[2,3:4]))
mfis_o18 <- melt(1- (mfis[6,3:4] / mfis[3,3:4]))



wf_mfi <- rbind(mfis_u10, mfis_o10,mfis_o18)
wf_mfi$Age_group <- c("COVID 1-10","COVID 1-10","COVID 11-18","COVID 11-18","COVID >18","COVID >18")
wf_mfi$Age_group <- factor(wf_mfi$Age_group, levels= c("COVID 1-10", "COVID 11-18", "COVID >18"))
wf_mfi$variable <- as.factor(wf_mfi$variable)
wf_mfi$variable <- str_remove(wf_mfi$variable, "_mfi_HLA-DR")
wf_mfi$variable <- str_replace_all(wf_mfi$variable, "_", " ")

wf_mfi$value <- (wf_mfi$value - (2*wf_mfi$value)) *100

wf_mfi$start <- 0
wf_mfi$id <- c(1,2,1,2,1,2)
wf_mfi$direction <- ifelse(wf_mfi$value > 0, "positive", "negative")



waterfall_covid_healthy <- wf_mfi %>% mutate(variable = forcats::fct_inorder(variable)) %>%
  ggplot() +
  facet_grid(.~ Age_group)+
  geom_rect(aes(x = variable, xmin = id -0.45, xmax = id+0.45,ymin = start, ymax= value,fill = direction
  ), color = "#252525", alpha = .7 , binwidth = 0.5) +
  geom_hline(yintercept = 0, alpha=0.8) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
  theme_classic() +
  scale_fill_manual(values = c(scales::muted("blue"),
                               scales::muted("red"))) + 
  labs(y = "% change of HLA-DR (MFI)
relative to healthy controls ") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust=1, size = 9),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text  = element_text(size=12, face="bold"),
        plot.caption = element_text(size = 10),
        axis.title.x = element_blank())

waterfall_covid_healthy


# DC MFI tile plot #####
dc_mfi_z <- as.data.frame(pall_clean_T0[,c("merge_ID","mDC_mfi_HLA-DR","pDC_mfi_HLA-DR" )])
rownames(dc_mfi_z) <- dc_mfi_z$merge_ID
# ex COV040 
rownames(dc_mfi_z)
dc_mfi_z <- dc_mfi_z[,-1]

#normalize into z-scores
for (i in 1:ncol(dc_mfi_z)) {
  m <- sum(dc_mfi_z[,i])/length(dc_mfi_z[,i])
  std <- sd(dc_mfi_z[,i])
  print(m)
  for (j in 1:nrow(dc_mfi_z)) {
    dc_mfi_z[j,i] <- (dc_mfi_z[j,i] - m) / std
  }
}

dc_mfi_z$merge_ID <- rownames(dc_mfi_z)
dc_mfi_z <- dc_mfi_z %>% left_join(pall_clean[,c("merge_ID","Status","Age_group")],by = "merge_ID")
colnames(dc_mfi_z) <- str_replace(colnames(dc_mfi_z), "-", "_")
colnames(dc_mfi_z)

mean_z_scores_dc <- aggregate(cbind(mDC_mfi_HLA_DR, pDC_mfi_HLA_DR) ~ Age_group + Status, 
                           data = dc_mfi_z, FUN = mean, na.rm = TRUE)
mean_z_scores_dc <- mean_z_scores_dc %>% gather(key = "parameter", value = "value", c(3:4))

mean_z_scores_dc$Status <- factor(mean_z_scores_dc$Status, levels=c("Non-COVID", "COVID", "Healthy"))
mean_z_scores_dc$parameter <- str_remove(mean_z_scores_dc$parameter, "_mfi_HLA_DR")
mean_z_scores_dc$parameter <- str_replace(mean_z_scores_dc$parameter, "_", " ")
mean_z_scores_dc$merge <-paste(mean_z_scores_dc$Age_group, mean_z_scores_dc$parameter, mean_z_scores_dc$Status)

pall_long$parameter <- str_replace(pall_long$parameter, "-", "_")
mfis_sign <- pall_long %>% filter(Status == "Healthy" | Status == "COVID" | Status =="Non-COVID") %>%
  filter(day == 0 & outlier_code < 5) %>%
  filter( Age_group != "<1") %>%
  filter(parameter == "mDC_mfi_HLA_DR"| parameter== "pDC_mfi_HLA_DR") %>% 
  group_by(parameter, Age_group)%>%
  wilcox_test(value  ~ Status)%>%
  adjust_pvalue(method = "BH") %>% 
  filter((group1 == "Healthy" & group2=="Non-COVID")|(group1 == "Healthy" & group2=="COVID")) %>%
  mutate(parameter = str_remove(parameter, "_mfi_HLA_DR")) %>% 
  mutate(parameter = str_replace(parameter, "_", " ")) %>% 
  mutate(merge = paste(Age_group, parameter, group2)) 

mfis_sign <- mfis_sign %>% 
  mutate(stars = cut(mfis_sign$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))) ## cut??

mean_z_scores_dc <- mean_z_scores_dc %>% left_join(mfis_sign[,9:13], by = "merge") 
mean_z_scores_dc$stars[is.na(mean_z_scores_dc$stars)] <- "" 
mean_z_scores_dc$p.adj[is.na(mean_z_scores_dc$p.adj)] <- 1
mean_z_scores_dc$FDR <- as.factor(ifelse(mean_z_scores_dc$p.adj < .1, 1,0))

class(mean_z_scores_dc$FDR)


ggplot(data = mean_z_scores_dc, aes(x=parameter, y=Status,fill=value)) + 
  facet_grid(.~Age_group,
             labeller = labeller(Age_group = p.labs.truncated)) +
  geom_tile(color="white", size = 1.2)+
  geom_tile(aes(color = FDR),fill = NA, size=1.7)+
  coord_equal() +
  geom_text(aes(label=stars), color="white", size=10) +
  scale_fill_gradient2(name = "HLA-DR
expression", 
                       breaks = c(.82,-0.87), labels = c("high", "low"),
                       midpoint = 0, low = scales::muted("blue"), mid = "gray96",
                       high = scales::muted("red"), space = "Lab", n.breaks = 2, 
                       guide = "colorbar")+
  scale_color_manual(name ="FDR", values=c("white","black"),
                     labels = c("0"="n.s","1"="FDR<0.1")) +
  #labs(x = "Cell populations",
  #    y= "Clinical measurments",
  #   fill = "Spearman's rho")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05,
                                   hjust=1),
        axis.text = element_text(size = 20),
        strip.background = element_blank(),
        strip.text = element_text(size =23, face="bold"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title =element_blank(),
        #legend.position = "right", 
        legend.text = element_text(size=16),
        legend.title = element_text(size=17))


ggsave("DC_Hla_Dr_tiles_1.pdf", height = 5, width = 7)

# line graph mono subpopulations hla-dr #####
#pdf("~/Ped-Covid Studie/FACS/Plots/line_graph_monocyte_subpopulations.pdf", onefile = T, paper="a4r",width=20, height=20)
ggplot(pall_long %>% filter(Age_group != "<1" &  outlier_code<5 & day == 0 & Status != "Non-COVID" &
                              (parameter == "non_class_mono_mfi_HLA-DR" | 
                                 parameter == "intermed_mono_mfi_HLA-DR" | 
                                 parameter ==  "class_mono_mfi_HLA-DR")), 
       aes(x =as.factor(parameter), y= value, color = Status, group = analysis.identifier.x)) +
  facet_grid(.~Age_group) +
  geom_point(aes(fill=Status), size=3, shape =21,alpha=.7, color="black", stroke =1.7) +
  geom_line(size = 1, alpha =.6) +
  scale_color_manual(values =c("#0381A1", "#BC5090")) +
  scale_fill_manual(values =c("#0381A1", "#BC5090")) +
  scale_x_discrete(labels=c("non_class_mono_mfi_HLA-DR" = "NC", "intermed_mono_mfi_HLA-DR" = "I",
                            "class_mono_mfi_HLA-DR" = "C"))+
  theme_classic2() +
  labs(y = "MFI HLA-DR") +
  theme( plot.title = element_text(size = 26, face = "bold"),
          axis.text = element_text( size = 20 ),
          #axis.text.x = element_text(angle = 45, vjust = 1.05, hjust=1),
          axis.title = element_text( size = 21),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text( size = 23, face="bold"),
          strip.background = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 18))

ggsave("mono_linegra.pdf", height = 6, width = 7)
  


# CD4 CD8 absolute  #####
pall_clean <- pall_clean %>% mutate(CD4_abs = ((CD4*CD3)/100)) %>%
  mutate(CD8_abs = ((CD8*CD3)/100)) %>% 
  mutate(Tregs_abs = ((Tregs*CD3)/100))

collapsed.1 <- pall_clean %>% mutate(Age_group = "All Patients")
collapsed.1 <- rbind(pall_clean, collapsed.1)

my_comparisons2 <- list(c("H","C"), c("C", "NC"))
collapsed.1$Age_group <- factor(collapsed.1$Age_group, levels =c("All Patients", "1-10","11-18",">18"))



collapsed_long <- collapsed.1 %>% gather(key = "parameter", value = "value", c(187:189))

p.labs <- c("CD4+", "CD8+")
names(p.labs) <- c("CD4_abs", "CD8_abs")


ggplot(collapsed_long %>% filter((parameter == "CD4_abs"| parameter == "CD8_abs") & day== 0), 
       aes(x = short_status, y = value, fill = short_status)) +
  facet_grid(.~parameter,
             labeller = labeller(parameter = p.labs)) +
  geom_boxplot(color = "black", outlier.shape = NA, alpha =0.4) +
  geom_jitter(aes(shape=Disease_course, fill = short_status), size=3.5, color = "black", stroke = 1.5,
              position = position_jitterdodge(jitter.width = 2.5, dodge.width = 0)) +
  scale_fill_manual(name = "COVID Status",
                    labels = c("Healthy","COVID","Non-COVID"),
                    values = c("#0381A1","#BC5090","#FFA600"))+
  scale_shape_manual(name = "Disease Severity ",
                     values = c(21:25)) +
  stat_compare_means(comparisons = my_comparisons2 ,aes(label=..p.adj..),
                     method = "wilcox.test", paired = FALSE, size = 8) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  #ggtitle("T-cells") +
  labs(y = "% / CD45+")+ 
  theme_linedraw() + 
  theme_linedraw_nogrid_facet  +
  theme( plot.title = element_text(size = 26, face = "bold"),
         axis.text = element_text( size = 20 ),
         #axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.title = element_text( size = 21),
         strip.text = element_text( size = 23, face="bold"),
         strip.background = element_blank(),
         legend.title = element_blank(),
         legend.text = element_text(size = 18),
         legend.position = "none",
         legend.key = element_rect(color= NA),
         legend.key.size = unit(1.0, "cm"),
         axis.ticks.x = element_blank(),
         legend.box="vertical", 
         legend.margin=margin())



ggsave("CD4CD8.pdf", height = 6, width = 4)





ggplot(collapsed_long %>% filter(parameter == "CD8_abs" & day== 0), 
       aes(x = Status, y = value, fill = Status)) +
  facet_grid(.~Age_group) +
  geom_boxplot(color = "black", outlier.shape = NA, alpha =0.4) +
  geom_jitter(aes(shape=Disease_course, fill = Status), size=3.5, color = "black", stroke = 1.5,
              position = position_jitterdodge(jitter.width = 2.5, dodge.width = 0)) +
  scale_fill_manual(name = "COVID Status",
                    labels = c("Healthy","COVID","Non-COVID"),
                    values = c("#0381A1","#BC5090","#FFA600"))+
  scale_shape_manual(name = "Disease Severity ",
                     values = c(21:25)) +
  stat_compare_means(comparisons = my_comparisons2 ,aes(label=..p.adj..),
                     method = "wilcox.test", paired = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  ggtitle("CD8+ T-cells") +
  labs(y = "% / Leukocytes")+ 
  theme_linedraw() + 
  theme_linedraw_nogrid_facet  +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.title = element_text( size = 14),
         strip.text = element_text( size = 16, face="bold"),
         strip.background = element_blank(),
         legend.title = element_text ( size = 13, face = "bold"),
         legend.text = element_text(size = 11),
         legend.position = "right",
         legend.key = element_rect(color= NA),
         legend.key.size = unit(1.0, "cm"),
         axis.ticks.x = element_blank(),
         legend.box="vertical", 
         legend.margin=margin())

ggsave("CD8_1.pdf", height = 6, width = 11)





ggplot(collapsed, 
       aes(x = Status, y = Tregs_abs, fill = Status)) +
  facet_grid(.~Age_group) + 
  geom_boxplot(color = "black", outlier.shape = NA, alpha =0.4) +
  geom_jitter(aes(shape=Disease_course, fill = Status), size=3.5, color = "black", stroke = 1.5,
              position = position_jitterdodge(jitter.width = 2.5, dodge.width = 0)) +
  scale_fill_manual(name = "COVID Status",
                    labels = c("Healthy","COVID","Non-COVID"),
                    values = c("#0381A1","#BC5090","#FFA600"))+
  scale_shape_manual(name = "Disease Severity ",
                     values = c(21:25)) +
  stat_compare_means(comparisons = my_comparisons2,aes(label=..p.adj..),
                     method = "wilcox.test", paired = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme_bw() +
  ggtitle("Tregs (abs.)") +
  theme( plot.title = element_text(size = 18, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.title = element_text( size = 14),
         strip.text = element_text( size = 16, face="bold"),
         strip.background = element_blank(),
         legend.title = element_text ( size = 13, face = "bold"),
         legend.text = element_text(size = 11),
         legend.position = "top",
         legend.key = element_rect(color= NA),
         legend.key.size = unit(1.0, "cm"),
         axis.ticks.x = element_blank(),
         legend.box="vertical", 
         legend.margin=margin())
dev.off()

# CD38 HLA-DR aktivierte Zellen  ####
pall_clean_T0 <- pall_clean %>% filter(day == 0 )
t_mfi_z <- as.data.frame(pall_clean_T0[-50,c(46,28:29,31:32)])
rownames(t_mfi_z) <- t_mfi_z$merge_ID
# ex COV040 
rownames(t_mfi_z)
t_mfi_z <- t_mfi_z[,-1]

#normalize into z-scores
for (i in 1:ncol(t_mfi_z)) {
  m <- sum(t_mfi_z[,i])/length(t_mfi_z[,i])
  std <- sd(t_mfi_z[,i])
  print(m)
  for (j in 1:nrow(t_mfi_z)) {
    t_mfi_z[j,i] <- (t_mfi_z[j,i] - m) / std
  }
}

t_mfi_z$merge_ID <- rownames(t_mfi_z)
t_mfi_z <- t_mfi_z %>% left_join(pall_clean[,c("Age_group", "Status","merge_ID")],by = "merge_ID")
colnames(t_mfi_z) <- str_replace(colnames(t_mfi_z), "-", "_")
colnames(t_mfi_z)

t_mean_z <- aggregate(cbind(activ_CD4m_mfi_HLA_DR,activ_CD8m_mfi_HLA_DR,
                            activ_CD4m_mfi_CD38, activ_CD8m_mfi_CD38) ~ Age_group + Status, 
                           data = t_mfi_z, FUN = mean, na.rm = TRUE)
t_mean_z <- t_mean_z %>% gather(key = "parameter", value = "value", c(3:6))

t_mean_z$Status <- factor(t_mean_z$Status, levels=c("Non-COVID", "COVID", "Healthy"))
t_mean_z$parameter <- str_replace(t_mean_z$parameter, "_mfi_", "/")
t_mean_z$parameter <- str_replace(t_mean_z$parameter, "_", " ")
t_mean_z$merge <-paste(t_mean_z$Age_group, t_mean_z$parameter, t_mean_z$Status)

pall_long$parameter <- str_replace(pall_long$parameter, "-", "_")
mfis_sign <- pall_long %>% filter(Status == "Healthy" | Status == "COVID" | Status =="Non-COVID") %>% 
  filter(day == 0 & outlier_code < 5) %>%
  filter( Age_group != "<1") %>%
  filter(parameter == "activ_CD4m_mfi_HLA_DR"|parameter == "activ_CD8m_mfi_HLA_DR"|
         parameter == "activ_CD4m_mfi_CD38"| parameter == "activ_CD8m_mfi_CD38") %>% 
  group_by(parameter, Age_group) %>%
  wilcox_test(value  ~ Status)%>%
  adjust_pvalue(method = "BH")  %>% 
  filter((group1 == "Healthy" & group2=="Non-COVID")|(group1 == "Healthy" & group2=="COVID")) %>%
  mutate(parameter = str_replace(parameter, "_mfi_", "/")) %>% 
  mutate(parameter = str_replace(parameter, "_", " ")) %>% 
  mutate(merge = paste(Age_group, parameter, group2)) 

mfis_sign <- mfis_sign %>% 
  mutate(stars = cut(mfis_sign$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))) ## cut??

t_mean_z <- t_mean_z %>% left_join(mfis_sign[,9:13], by = "merge") 
t_mean_z$stars[is.na(t_mean_z$stars)] <- "" 
t_mean_z$p.adj[is.na(t_mean_z$p.adj)] <- 1
t_mean_z$FDR <- as.factor(ifelse(t_mean_z$p.adj < .1, 1,0))
t_mean_z$parameter <- factor(t_mean_z$parameter, levels=c("activ CD4m/CD38", "activ CD8m/CD38",
                                                          "activ CD4m/HLA_DR","activ CD8m/HLA_DR"))
class(t_mean_z$FDR)


ggplot(data = t_mean_z, aes(x=parameter, y=Status, 
                                 fill=value)) + 
  facet_grid(.~Age_group,
             labeller = labeller(Age_group = p.labs.truncated)) +
  geom_tile(color = "white", size = 1.2)+
  geom_tile(aes(color = FDR),fill = NA, size=1.7)+
  coord_equal() +
  geom_text(aes(label=stars), color="white", size=8) +
  scale_fill_gradient2(name = "Expression
level", #labels = c("1" = "high", "-1" = "low"),
                       breaks = c(1.78,-0.9), labels = c("high", "low"),
                       midpoint = 0, low = scales::muted("blue"), mid = "gray96",
                       high = scales::muted("red"), space = "Lab", n.breaks = 2)+
  scale_color_manual(name ="FDR", values=c("white","black"),
                     labels = c("0"="n.s","1"="FDR<0.1")) +
  #labs(x = "Cell populations",
  #    y= "Clinical measurments",
  #   fill = "Spearman's rho")+
  theme_classic() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, 
         #                          hjust=1, size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        strip.background = element_blank(),
        strip.text = element_text(size =23, face="bold"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title =element_blank(),
        legend.position = "right", 
        legend.text = element_text(size=16),
        legend.title = element_text(size=17))

ggsave("CD38HLADR_Tcells_excl_labels.pdf", height = 6, width = 10)


# B cells absolute ######


pall_clean <- pall_clean %>% mutate(naive_CD19_abs = ((CD19*naive_CD19)/100)) %>%
  mutate(NCS_CD19_abs = ((CD19*non_switch_CD19m)/100)) %>% 
  mutate(CS_CD19_abs = ((CD19*class_switch_CD19m)/100)) %>% 
  mutate(trans_abs =((transitional*class_switch_CD19m)/100) ) %>% 
  mutate(pbs_abs =((plasmablasts*CD19)/100) )

collapsed <- pall_clean %>% mutate(Age_group = "All Patients")
collapsed <- rbind(pall_clean, collapsed)

my_comparisons2 <- list(c("Healthy","COVID"), c("COVID", "Non-COVID"))
collapsed$Age_group <- factor(collapsed$Age_group, levels =c("All Patients", "1-10","11-18",">18"))




ggplot(collapsed %>% filter(day== 0 & NCS_CD19_abs < 4), 
       aes(x = Status, y = NCS_CD19_abs, fill = Status)) +
  facet_grid(.~Age_group) + 
  geom_boxplot(color = "black", outlier.shape = NA, alpha =0.4) +
  geom_jitter(aes(shape=Disease_course, fill = Status), size=3.5, color = "black", stroke = 1.5,
              position = position_jitterdodge(jitter.width = 2.5, dodge.width = 0)) +
  scale_fill_manual(name = "COVID Status",
                    labels = c("Healthy","COVID","Non-COVID"),
                    values = c("#0381A1","#BC5090","#FFA600"))+
  scale_shape_manual(name = "Disease Severity ",
                     values = c(21:25)) +
  stat_compare_means(comparisons = my_comparisons2,aes(label=..p.adj..),
                     method = "wilcox.test", paired = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  ggtitle("Non-switch memory B cells") +
  labs(y = "% / CD45+")+ 
  theme_linedraw() + 
  theme_linedraw_nogrid_facet  +
  theme( plot.title = element_text(size = 22, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.title = element_text( size = 20),
         strip.text = element_text( size = 16, face="bold"),
         strip.background = element_blank(),
         legend.title = element_text ( size = 13, face = "bold"),
         legend.text = element_text(size = 11),
         legend.position = "none",
         legend.key = element_rect(color= NA),
         legend.key.size = unit(1.0, "cm"),
         axis.ticks.x = element_blank(),
         legend.box="vertical", 
         legend.margin=margin())

ggsave("NCS_Bcells.pdf", height = 6, width = 7)


ggplot(collapsed %>% filter(day== 0 ), 
       aes(x = Status, y = CS_CD19_abs, fill = Status)) +
  facet_grid(.~Age_group) + 
  geom_boxplot(color = "black", outlier.shape = NA, alpha =0.4) +
  geom_jitter(aes(shape=Disease_course, fill = Status), size=3.5, color = "black", stroke = 1.5,
              position = position_jitterdodge(jitter.width = 2.5, dodge.width = 0)) +
  scale_fill_manual(name = "COVID Status",
                    labels = c("Healthy","COVID","Non-COVID"),
                    values = c("#0381A1","#BC5090","#FFA600"))+
  scale_shape_manual(name = "Disease Severity ",
                     values = c(21:25)) +
  stat_compare_means(comparisons = my_comparisons2,aes(label=..p.adj..),
                     method = "wilcox.test", paired = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  ggtitle("Class-switch memory B cells") +
  labs(y = "% / CD45+")+ 
  theme_linedraw() + 
  theme_linedraw_nogrid_facet  +
  theme( plot.title = element_text(size = 22, face = "bold"),
         axis.text = element_text( size = 16 ),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.title = element_text( size = 20),
         strip.text = element_text( size = 16, face="bold"),
         strip.background = element_blank(),
         legend.title = element_text ( size = 13, face = "bold"),
         legend.text = element_text(size = 11),
         legend.position = "none",
         legend.key = element_rect(color= NA),
         legend.key.size = unit(1.0, "cm"),
         axis.ticks.x = element_blank(),
         legend.box="vertical", 
         legend.margin=margin())

ggsave("CS_Bcells.pdf", height = 6, width = 7)


# tiles Monos MFI CD11c #######
mono_mfi_z <- as.data.frame(pall_clean_T0[,c("merge_ID","classs_mono_cd11c_mean","intermed_mono_cd11c_mean",
                                          "non_class_mono_cd11c_mean")])
rownames(mono_mfi_z) <- mono_mfi_z$merge_ID
# ex COV040 
rownames(mono_mfi_z)
mono_mfi_z <- mono_mfi_z[-c(31,46,51,54,55),-1]

#normalize into z-scores
for (i in 1:ncol(mono_mfi_z)) {
  m <- sum(mono_mfi_z[,i])/length(mono_mfi_z[,i])
  std <- sd(mono_mfi_z[,i])
  print(m)
  for (j in 1:nrow(mono_mfi_z)) {
    mono_mfi_z[j,i] <- (mono_mfi_z[j,i] - m) / std
  }
}

mono_mfi_z$merge_ID <- rownames(mono_mfi_z)
mono_mfi_z <- mono_mfi_z %>% left_join(pall_clean[,c("merge_ID","Status","Age_group")],by = "merge_ID")
colnames(mono_mfi_z)

mean_z_scores <- aggregate(cbind(classs_mono_cd11c_mean,intermed_mono_cd11c_mean,
                                 non_class_mono_cd11c_mean) ~ Age_group + Status, 
                           data = mono_mfi_z, FUN = mean, na.rm = TRUE)
mean_z_scores <- mean_z_scores %>% gather(key = "parameter", value = "value", c(3:5))
mean_z_scores$Status <- factor(mean_z_scores$Status, levels=c("Non-COVID", "COVID", "Healthy"))

mean_z_scores$merge <-paste(mean_z_scores$Age_group, mean_z_scores$parameter, mean_z_scores$Status)

unique(pall_long$parameter)

pall_long$parameter <- str_replace(pall_long$parameter, "-", "_")
mfis_sign <- pall_long %>% filter(Status == "Healthy" | Status == "COVID" | Status =="Non-COVID") %>%
  filter(day == 0 & outlier_code < 5) %>%
  filter( Age_group != "<1") %>%
  filter(parameter == "classs_mono_cd11c_mean"| parameter== "intermed_mono_cd11c_mean"|
           parameter =="non_class_mono_cd11c_mean") %>% 
  group_by(parameter, Age_group) %>%
  wilcox_test(value  ~ Status)%>%
  adjust_pvalue(method = "BH")  %>% 
    filter((group1 == "Healthy" & group2=="Non-COVID")|(group1 == "Healthy" & group2=="COVID")) %>% 
  mutate(merge = paste(Age_group, parameter, group2))

mfis_sign <- mfis_sign %>% 
  mutate(stars = cut(mfis_sign$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))) ## cut??

mean_z_scores <- mean_z_scores %>% left_join(mfis_sign[,9:13], by = "merge") 
mean_z_scores$stars[is.na(mean_z_scores$stars)] <- "" 
mean_z_scores$p.adj[is.na(mean_z_scores$p.adj)] <- 1
mean_z_scores$FDR <- as.factor(ifelse(mean_z_scores$p.adj < .1, 1,0))

mean_z_scores$parameter[mean_z_scores$parameter == "classs_mono_cd11c_mean"] <- "Classical"
mean_z_scores$parameter[mean_z_scores$parameter == "intermed_mono_cd11c_mean"] <- "Intermediate"
mean_z_scores$parameter[mean_z_scores$parameter == "non_class_mono_cd11c_mean"] <- "Non-Classical"


class(mean_z_scores$FDR)


ggplot(data = mean_z_scores, aes(x=parameter, y=Status, 
                                 fill=value)) + 
  facet_grid(.~Age_group,
             labeller = labeller(Age_group = p.labs.truncated)) +
  geom_tile(color = "white", size = 1.2)+
  geom_tile(aes(color = FDR),fill = NA, size=1.7)+
  coord_equal() +
  geom_text(aes(label=stars), color="white", size=10) +
  scale_fill_gradient2(name = "CD11c
expression", #labels = c("1" = "high", "-1" = "low"),
                       breaks = c(1.4,-0.6), labels = c("high", "low"),
                       midpoint = 0, low = scales::muted("blue"), mid = "gray96",
                       high = scales::muted("red"), space = "Lab", n.breaks = 2)+
  scale_color_manual(name ="FDR", values=c("white","black"),
                     labels = c("0"="n.s","1"="FDR<0.1")) +
  #labs(x = "Cell populations",
  #    y= "Clinical measurments",
  #   fill = "Spearman's rho")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05,
                                   hjust=1),
        axis.text = element_text(size = 20),
        strip.background = element_blank(),
        strip.text = element_text(size =23, face="bold"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title =element_blank(),
        legend.position = "right", 
        legend.text = element_text(size=16),
        legend.title = element_text(size=17))


ggsave("Mono_CD11c_tiles_large.pdf", height = 6, width = 10)

# tiles mDCs MFI CD11c #########

dc_mfi_z <- as.data.frame(pall_clean_T0[,c("merge_ID", "mdc_cd11c_mean")])
rownames(dc_mfi_z) <- dc_mfi_z$merge_ID
# ex COV040 
rownames(dc_mfi_z)
dc_mfi_z <- dc_mfi_z[-c(46,51,54,55),]

#normalize into z-scores
  m <- sum(dc_mfi_z[,2])/length(dc_mfi_z[,2])
  std <- sd(dc_mfi_z[,2])
  print(m)
  for (j in 1:length(dc_mfi_z[,2])) {
    dc_mfi_z[j,2] <- (dc_mfi_z[j,2] - m) / std}
  

dc_mfi_z <- dc_mfi_z %>% left_join(pall_clean[,c("merge_ID", "Status", "Age_group")],by = "merge_ID")
colnames(dc_mfi_z) <- str_replace(colnames(dc_mfi_z), "-", "_")
colnames(dc_mfi_z)

mean_z_scores_dc <- aggregate(mdc_cd11c_mean ~ Age_group + Status, 
                              data = dc_mfi_z, FUN = mean, na.rm = TRUE)
mean_z_scores_dc <- mean_z_scores_dc %>% gather(key = "parameter", value = "value", c(3))

mean_z_scores_dc$Status <- factor(mean_z_scores_dc$Status, levels=c("Non-COVID", "COVID", "Healthy"))

mean_z_scores_dc$merge <-paste(mean_z_scores_dc$Age_group, mean_z_scores_dc$parameter, mean_z_scores_dc$Status)

mfis_sign <- pall_long %>% filter(Status == "Healthy" | Status == "COVID" | Status =="Non-COVID") %>%
  filter(day == 0 & outlier_code < 5) %>%
  filter( Age_group != "<1") %>%
  filter(parameter == "mdc_cd11c_mean") %>% 
  group_by(parameter, Age_group)%>%
  wilcox_test(value  ~ Status)%>%
  adjust_pvalue(method = "BH") %>% 
  filter((group1 == "Healthy" & group2=="Non-COVID")|(group1 == "Healthy" & group2=="COVID")) %>%
  mutate(merge = paste(Age_group, parameter, group2)) 

mfis_sign <- mfis_sign %>% 
  mutate(stars = cut(mfis_sign$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))) ## cut??

mean_z_scores_dc <- mean_z_scores_dc %>% left_join(mfis_sign[,9:13], by = "merge") 
mean_z_scores_dc$stars[is.na(mean_z_scores_dc$stars)] <- "" 
mean_z_scores_dc$p.adj[is.na(mean_z_scores_dc$p.adj)] <- 1
mean_z_scores_dc$FDR <- as.factor(ifelse(mean_z_scores_dc$p.adj < .1, 1,0))

class(mean_z_scores_dc$FDR)

mean_z_scores_dc$parameter <- "mDC"


ggplot(data = mean_z_scores_dc, aes(x=parameter, y=Status,fill=value)) + 
  facet_grid(.~Age_group,
             labeller = labeller(Age_group = p.labs.truncated)) +
  geom_tile(color="white", size = 1.2)+
  geom_tile(aes(color = FDR),fill = NA, size=1.7)+
  coord_equal() +
  geom_text(aes(label=stars), color="white", size=10) +
  scale_fill_gradient2(name = "CD11c
expression", 
                       breaks = c(.92,-0.72), labels = c("high", "low"),
                       midpoint = 0, low = scales::muted("blue"), mid = "white",
                       high = scales::muted("red"), space = "Lab", n.breaks = 2)+
  scale_color_manual(name ="FDR", values=c("white","black"),
                     labels = c("0"="n.s","1"="FDR<0.1")) +
  #labs(x = "Cell populations",
  #    y= "Clinical measurments",
  #   fill = "Spearman's rho")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05,
                                   hjust=1),
        axis.text = element_text(size = 20),
        strip.background = element_blank(),
        strip.text = element_text(size =19, face="bold"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title =element_blank(),
        legend.position = "right", 
        legend.text = element_text(size=16),
        legend.title = element_text(size=17))

ggsave("mDC_CD11c_tiles_2.pdf", height = 4, width = 5)
