#libraries
#####
library(purrr) # v. 0.3.4
library(ggplot2) # v. 3.3.2
library(patchwork) # v. 1.0.0
library(broom)
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
library(PCAtools)
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
#####

setwd("~/Ped-Covid Studie")
load("~/Ped-Covid Studie/FACS/scripts/all_panels_common_WS.RData")

complete_panels <- Reduce(intersect, list(panel1_cov$merge_ID,panel2_cov$merge_ID,P3_cov_complete$merge_ID))

p1_common <- panel1_cov[panel1_cov$merge_ID %in% complete_panels,]
p2_common <- panel2_cov[panel2_cov$merge_ID %in% complete_panels,]
p3_common <- P3_cov_complete[P3_cov_complete$merge_ID %in% complete_panels,]

#merge all panels together, exclude outliers , make second timepoint of COV068 day 0 because PCR tst was positive for the first time only before that
pall <- p1_common[,c(6:10,12,13,15,17:19,21,22,25:41,52,53,55,57,58,60,61,63,64,66,67,76,77,79,80,172)] %>% 
  left_join(p2_common[,c(6:11,13:15,98)], by = "merge_ID") %>% 
  left_join(p3_common[,c(4:19,24:29,33,116)], by = "merge_ID") %>% 
  left_join(p1_common[,90:181], by = "merge_ID") %>%
  mutate(day = replace(day, merge_ID == 	"COV068_1", -5)) %>%
  mutate(day = replace(day, merge_ID == 	"COV068_2", 0)) 

rm(list = c("p1_common", "p2_common", "p3_common", "P3_cov_complete", "panel1_cov", "panel2_cov"))

# add columns PCR symptom begin and COVID, outlier, sex dummy etc
pall$days_symptom_start <- pall$sampling_date - pall$date_of_1st_symptoms
pall$days_symptom_start <- as.numeric(str_remove(pall$days_symptom_start, " days"))

pall$days_pcr_sampling <- pall$sampling_date - pall$date_of_1st_positive_PCR
pall$days_pcr_sampling <- as.numeric(str_remove(pall$days_pcr_sampling, " days"))

pall$OC5 <- ifelse(pall$outlier_code == 5, 1, 0)
pall$sex01 <- pall$`Sex_(f=0;_m=1)`
pall$covid <- ifelse(pall$Status == "COVID", 1, 0)
pall$late_sampling <- as.factor(ifelse(pall$days_pcr_sampling > 7, 1, 0))

#exclude MFI columns 
pattern_mfi <- or("MFI", "mfi")
pall_exMFI <- pall %>% select(all_of(colnames(pall)[!str_detect(colnames(pall), pattern_mfi)]))

#two data sets: COVID & Healhty+COVID
pall_cov <- pall_exMFI %>% filter(Status == "COVID"  & day == 0)
pall_cov_h <- pall_exMFI %>% filter((Status == "COVID" | Status == "Healthy") & day == 0 )


# compile function for modelling & evaluating 
model_fit_test <- function(panel, ind_variable, filename, print_resp = FALSE){
  ttest_fun = function(response) {
    form = paste(response, ind_variable)
    lm(as.formula(form), data = panel )
  }
  
  ###summary(ttest_fun("naive_CD8"))
  vars <- names(panel)[c(1:27,29:53)]
  
  models <- vars %>% set_names() %>% map(ttest_fun)  # map: Apply a function to each element of a list or atomic vector # set_nameS(): sets names for vecotr x, name vector must be same lenght as x 
  
  #create residual plots in a singel picture using patchwork 
  resid_plots <- function(model, modelname) {
    output = augment(model) # augment: adds information: fitted values, residuals & standard error 
    
    res.v.fit = ggplot(output, aes(x = .fitted, y = .resid) ) +
      geom_point() +
      theme_bw(base_size = 16)
    
    res.box = ggplot(output, aes(x = "", y = .resid) ) +
      geom_boxplot() +
      theme_bw(base_size = 16) +
      labs(x = NULL)
    
    res.v.fit + res.box +
      plot_annotation(title = paste("Residuals plots for", modelname) )
  }
  #resid_plots(model=models[[1]], modelname = names(models)[1])
  residplots <- imap(models, resid_plots)
  if (print_resp == TRUE) {
    pdf(paste("~/Päd-Covid Studie/FACS/Plots/final/", filename, "_residualplots.pdf", sep =""))
    print(residplots)
    dev.off()
  }
  #
  res_anova <- map_dfr(models, tidy, conf.int = FALSE, .id = "variable")
  res_anova_sign <- res_anova %>% filter(p.value < .05)
  # using all coefficients & only sign. ones
  write_xlsx(res_anova_sign, (paste(filename, "_sign_coefficients.xlsx", sep ="")))
  write_xlsx(res_anova, (paste(filename, "_all_coefficients.xlsx", sep ="")))
  print(res_anova)
}

model_fit_test(panel = pall_cov, ind_variable ="~ Age + WHO_on_date + OC5 + sex01  + days_pcr_sampling", filename="p1_covid_1" ) %>% view()

#####
## no exclusion critieria 
# COVID
model_fit_test(panel = pall_cov, ind_variable ="~ Age + WHO_on_date + OC5 + sex01  + days_pcr_sampling", filename="p1_covid_1" )
# COVID & Healthy 
model_fit_test(panel = pall_cov_h, ind_variable ="~ Age + WHO_on_date + OC5 + sex01 + covid", filename="p1_covid_healthy_1" )

## ex < 1 y.o.
pall_cov_exU1 <- pall_cov %>% filter(Age > 0)
pall_cov_h_exU1 <- pall_cov_h %>% filter(Age > 0)
# COVID
model_fit_test(panel = pall_cov_exU1, ind_variable ="~ Age + WHO_on_date + OC5 + sex01  + days_pcr_sampling", filename="p1_covid_2" )
# COVID & Healthy 
model_fit_test(panel = pall_cov_h_exU1, ind_variable ="~ Age + WHO_on_date + OC5 + sex01 + covid", filename="p1_covid_healthy_2" )

## ex > 18
pall_cov_exO18 <- pall_cov %>% filter(Age <= 18)
pall_cov_h_exO18 <- pall_cov_h %>% filter(Age <= 18)
# COVID
model_fit_test(panel = pall_cov_exO18, ind_variable ="~ Age + WHO_on_date + OC5 + sex01  + days_pcr_sampling", filename="p1_covid_3" )
# COVID & Healthy 
model_fit_test(panel = pall_cov_h_exO18, ind_variable ="~ Age + WHO_on_date + OC5 + sex01 + covid", filename="p1_covid_healthy_3" )

## ex <1 & >18
pall_cov_exage <- pall_cov %>% filter(Age <= 18 & Age >0)
pall_cov_h_exage <- pall_cov_h %>% filter(Age <= 18 & Age >0)
# COVID
model_fit_test(panel = pall_cov_exage, ind_variable ="~ Age + WHO_on_date + OC5 + sex01  + days_pcr_sampling", filename="p1_covid_4" )
# COVID & Healthy 
model_fit_test(panel = pall_cov_h_exage, ind_variable ="~ Age + WHO_on_date + OC5 + sex01 + covid", filename="p1_covid_healthy_4" )

## ex < 1 & > 18 & OC5
pall_cov_clean <- pall_cov %>% filter(Age <= 18 & Age >0 & outlier_code != 5)
pall_cov_h_clean <- pall_cov_h %>% filter(Age <= 18 & Age >0 & outlier_code != 5)
# COVID
model_fit_test(panel = pall_cov_clean, ind_variable ="~ Age + WHO_on_date + OC5 + sex01  + days_pcr_sampling", filename="p1_covid_5" )
# COVID & Healthy 
model_fit_test(panel = pall_cov_h_clean, ind_variable ="~ Age + WHO_on_date + OC5 + sex01 + covid", filename="p1_covid_healthy_5" )

#####

pall_cov$update_status <- "complete"
pall_cov_exage$update_status <- "exAge"
pall_cov_clean$update_status <- "exAgeOC5"

pall_comp <- rbind(pall_cov, pall_cov_exage, pall_cov_clean)
pall_comp <- pall_comp %>%  gather(key = "parameter", value = "value", c(1:27,29:53))
pall_comp$update_status <- as.factor(pall_comp$update_status)

number_of_facets <- ceiling(length(levels(as.factor(pall_comp$parameter))))
parameters <- levels(as.factor(pall_comp$parameter))

qcplot <- ggplot(pall_comp, aes(x = parameter, y = value , fill = factor(update_status), color = factor(update_status))) + 
  geom_boxplot() + 
  theme_bw() +
  labs(title ="QC plot", 
       y = "%") +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(), 
        legend.title = element_blank())

pdf("~/Päd-Covid Studie/FACS/Plots/final/QC_2/20210504_panels_report_boxplots_comparison.pdf", onefile = T, paper="a4r",width=20, height=10)
print(qcplot)
for (i in 1:number_of_facets){
  print(
    ggplot(pall_comp %>%filter (parameter == parameters[i]), 
           aes(x=update_status, y=value, fill = update_status)) + 
      facet_grid(.~ Age_group) +
      geom_boxplot() + 
      geom_point(aes(shape = late_sampling, color = late_sampling), size = 3.5)+
      theme_classic() +
      geom_text(aes(label=merge_ID),hjust=1, vjust=0, size = 3) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
      ggtitle(parameters[i]) +
      theme( plot.title = element_text(size = 18, face = "bold"),
             axis.text = element_text( size = 16 ),
             axis.title = element_text( size = 16, face = "bold" ),
             legend.title = element_text(size = 14),
             legend.text = element_text(size = 12))
  )
  print(i)
}
dev.off()



pall_clean <- pall %>% filter(day == 0 & Age >0 & Age <= 18 & outlier_code < 5)
pall_clean <- pall_clean %>% select(all_of(colnames(pall_clean)[!str_detect(colnames(pall_clean), pattern_mfi)]))
pall_clean_long <- pall_clean %>% gather(key="parameter", value = "value", c(1:27,29:53))
pall_clean_long$Status <- factor(pall_clean_long$Status , levels=c("Healthy", "Non-COVID", "COVID"))


parameters <- levels(as.factor(pall_clean_long$parameter))
my_comparisons <- list( c("COVID", "Non-COVID"), c("COVID", "Healthy"), c("Non-COVID", "Healthy")
                       # ,c("COVID", "MISC"),c("Non-COVID", "MISC"),c("Healthy", "MISC")
                       )

pall_clean_long 


pdf("~/Päd-Covid Studie/FACS/Plots/final/QC_2/20210504_panels_clean_updated_exU1O25.pdf", onefile = T, paper="a4r",width=20, height=10)
for (i in 1:number_of_facets){
  print(
    ggplot(pall_clean_long %>% filter (parameter == parameters[i]), 
           aes(x=Status, y=value, fill = Status)) + 
      facet_grid(.~ Age_group) +
      geom_boxplot() + 
      geom_point(aes(shape = late_sampling, color = late_sampling), size = 3.5)+
      theme_classic() +
      #geom_text(aes(label=merge_ID),hjust=1, vjust=0, size = 3) +
      scale_fill_manual(labels = c("Healthy","Non-Covid","Covid"),
      values = c("#0381A1","#FFA600","#BC5090"))+
      stat_compare_means(comparisons = my_comparisons,aes(label = ..p.signif..),
                        method = "wilcox.test", paired = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
      ggtitle(parameters[i]) +
      theme( plot.title = element_text(size = 18, face = "bold"),
             axis.text = element_text( size = 16 ),
             axis.title = element_text( size = 16, face = "bold" ),
             legend.title = element_text(size = 14),
             legend.text = element_text(size = 12))
  )
  print(i)
}
dev.off()

aggregate(pall_clean$WHO_on_date, by = list(pall_clean$sex01, pall_clean$Status), mean)


lm_eqn <- function(df, y, formula, x){
  m <- lm(as.formula(formula), df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

ratioplot1
ggplot(pall_clean %>% filter(Status == "COVID") , aes(x=WHO_on_date)) + 
  geom_point(aes(y=CD4, color = "red")) + 
  geom_smooth(aes(y=CD4, color = "red"),method = "glm", formula = y~x) + 
  geom_text(x = 3, y = 70,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), CD4,"CD4 ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5)+
  geom_point(aes(y= CD8, color= "blue")) + 
  geom_smooth(aes(y=CD8, color= "blue"),method = "glm", formula = y~x) +
  geom_text(x = 3, y = 18,label = lm_eqn(pall_clean %>% filter(Status == "COVID"), CD8,"CD8 ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5)+
  labs(title = "CD4 & CD8 / WHO_on_Date",
       subtitle = "data = only COVID, Age >0 & <25, no outliers ",
       y = "%  / CD3+",
       caption = "CD4 = blue; CD8 = red") +
  theme_classic()


ggplot(pall %>% filter(Status == "COVID" & day == 0) , aes(x=WHO_on_date)) + 
  geom_point(aes(y=CD4, color = "red")) + 
  geom_smooth(aes(y=CD4, color = "red"),method = "glm", formula = y~x) + 
  geom_text(x = 3, y = 70,label = lm_eqn(pall %>% filter(Status == "COVID" & day == 0), 
                                         CD4,"CD4 ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5)+
  geom_point(aes(y= CD8, color= "blue")) + 
  geom_smooth(aes(y=CD8, color= "blue"),method = "glm", formula = y~x) +
  geom_text(x = 3, y = 18,label = lm_eqn(pall %>% filter(Status == "COVID" & day == 0), 
                                         CD8,"CD8 ~ WHO_on_date", WHO_on_date), parse = TRUE, size = 5)+
  labs(title = "CD4 & CD8 / WHO_on_Date",
       subtitle = "data = only COVID, Age >0 & <25, no outliers ",
       y = "%  / CD3+",
       caption = "CD4 = blue; CD8 = red") +
  theme_classic()





pall$Age_group[pall$Age >18] <- ">18"
pall$Age_group[pall$Age <1] <- "<1"

pall$Age_group <- factor(pall$Age_group , levels=c("<1", "<=10", ">10", ">18"))

pall_clean <- pall %>% filter(day == 0 & outlier_code < 5)

pall_clean_long <- pall_clean %>% gather(key="parameter", value = "value", c(1:45,47:78))
pall_clean_long$Status <- factor(pall_clean_long$Status , levels=c("Healthy", "Non-COVID", "COVID"))
pall_clean_long$Disease_course <- ifelse(pall_clean_long$max_WHO_classification >= 3, "moderate", "mild")
pall_clean_long$Disease_course[pall_clean_long$max_WHO_classification > 4] <- "severe"
pall_clean_long$Disease_course <- as.factor(pall_clean_long$Disease_course)

parameters <- levels(as.factor(pall_clean_long$parameter))
my_comparisons <- list( c("COVID", "Non-COVID"), c("COVID", "Healthy"), c("Non-COVID", "Healthy")
                        # ,c("COVID", "MISC"),c("Non-COVID", "MISC"),c("Healthy", "MISC")
)
number_of_facets <- ceiling(length(levels(as.factor(pall_clean_long$parameter))))

pdf("~/Päd-Covid Studie/FACS/Plots/final/20210507_FACS2.pdf", onefile = T, paper="a4r",width=20, height=10)
for (i in 1:number_of_facets){
  print(
    ggplot(pall_clean_long %>% filter (parameter == parameters[i]), 
           aes(x=Status, y=value, fill = Status)) + 
      facet_grid(.~ Age_group) +
      geom_boxplot() +
      #geom_jitter(aes(shape=Disease_course), alpha = 0.6, size=2.5,
      #           position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.4)) +
      ggbeeswarm::geom_beeswarm(aes(shape=Disease_course)) +
      ggrepel::geom_text_repel(aes(label=merge_ID), size = 2) +
      #geom_point(aes(shape = late_sampling, color = late_sampling), size = 3.5)+
      theme_classic() +
      #geom_text(aes(label=merge_ID),hjust=1, vjust=0, size = 3) +
      scale_fill_manual(labels = c("Healthy","Non-Covid","Covid"),
                        values = c("#0381A1","#FFA600","#BC5090"))+
      stat_compare_means(comparisons = my_comparisons,aes(label = ..p.signif..),
                         method = "wilcox.test", paired = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
      ggtitle(parameters[i]) +
      theme( plot.title = element_text(size = 18, face = "bold"),
             axis.text = element_text( size = 16 ),
             axis.title = element_text( size = 16, face = "bold" ),
             legend.title = element_text(size = 14),
             legend.text = element_text(size = 12))
    )
  print(i)
}
dev.off()


# lm for analysis  ####
pall_analysis <- pall %>% filter(Age_group != "<1" & Age_group != ">18" & outlier_code <5 & day== 0 ) %>% 
  select(all_of(colnames(pall)[!str_detect(colnames(pall), pattern_mfi)])) %>% 
  filter(Status == "Healthy" | Status == "COVID")


model_fit_test(panel = pall_analysis, ind_variable ="~Age + covid", filename="pall_age_covid_nointeraction" )


