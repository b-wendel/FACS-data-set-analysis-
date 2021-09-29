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


setwd("~/Ped-Covid Studie")
### COV046 - CD56 staining failed - exclude all NK and NKT cells 
### COV020 + COV040 - CD14 staining failed - exclude all cells downstream of Monocytes, pDC's and mDC's
#load dataset with correct col class

P3_all_samples <- read_delim("~/Ped-Covid Studie/FACS/panel3_freeze123_20210510.csv", 
                             ";", escape_double = FALSE, 
                             col_types = cols(.default = "n", sample_ID = col_character()), 
                             locale = locale(decimal_mark = ","), 
                             trim_ws = TRUE)

meta <- read_delim("~/Ped-Covid Studie/metadaten/patient_metadata_20210511.csv",
                   ";", escape_double = FALSE, col_types = cols(DOB = col_date(format = "%d.%m.%Y"),
                                                                Status = col_factor(levels = c("COVID","Non-COVID", "Healthy","MISC")), 
                                                                date_inclusion = col_date(format = "%d.%m.%Y"), 
                                                                date_of_1st_fever = col_date(format = "%d.%m.%Y"),
                                                                date_of_1st_negative_PCR = col_date(format = "%d.%m.%Y"),
                                                                date_of_1st_positive_AB_test = col_date(format = "%d.%m.%Y"),
                                                                date_of_1st_positive_PCR = col_date(format = "%d.%m.%Y"),
                                                                date_of_1st_symptoms = col_date(format = "%d.%m.%Y"),
                                                                date_of_ICU_admission = col_date(format = "%d.%m.%Y"),
                                                                date_of_contact = col_date(format = "%d.%m.%Y"),
                                                                date_of_hospital_arrival = col_date(format = "%d.%m.%Y"),
                                                                date_of_hospital_discharge = col_date(format = "%d.%m.%Y")),
                   trim_ws = TRUE)
meta <- meta[,c(1:27, 48:65,66,70,74,78,82,86,90,94,98,102,111:115,124:133,172:176)]

#create a column containign age groups 
meta$Age_group <- ifelse(meta$Age<=10, "1-10", "11-18")
meta$Age_group[meta$Age<1] <- "<1"
meta$Age_group[meta$Age>18] <- ">18"
meta$Age_group <- factor(meta$Age_group, levels = c("<1","1-10","11-18",">18"))
#meta$Age_group[meta$Age>18] <- ">18"
meta$sex <- ifelse(meta$`Sex_(f=0;_m=1)`== 0,"f", "m")

colnames(meta)[13] <- "outlier_code"


P3_all_samples %>% count(sample_ID) %>% 
  filter(n > 1)


duplicated(P3_all_samples)
sum(duplicated(P3_all_samples))

#add analysis.identifier column, remove Specimen_001_, name it P3_all_samples_clean

#patterns to identify COV and adult HD

pattern_1 <- or("COV", "HD")
str_view(P3_all_samples$sample_ID, pattern_1)
str_match(P3_all_samples$sample_ID, pattern_1)

#create column "subgroup": COV, HD

P3_all_samples_clean <- P3_all_samples %>%
  mutate(analysis.identifier = str_remove(sample_ID, "Specimen_001_")) %>%
  mutate(subgroup = str_extract(analysis.identifier, pattern_1))

#get rid of IBD/SY samples :-)
pattern_IBD_SY <- or("IBD", "SY", "ID35","HH","WholeBlood", "Mother", "PatTGO-Stuttgart")
str_view(P3_all_samples$sample_ID, pattern_IBD_SY)

P3_all_samples_clean <- P3_all_samples_clean %>%
  filter(!str_detect(sample_ID, pattern_IBD_SY))


P3_all_samples_clean <- P3_all_samples_clean %>%
  mutate(subgroup = replace(subgroup, which(is.na(subgroup)), "time course"))

#simplify analysis identifier to "Cov-XX"

pattern_2 <- "_" %R% DGT %R% DGT %R% DGT %R% ".fcs"
str_view(P3_all_samples_clean$analysis.identifier, pattern_2)
P3_all_samples_clean <- P3_all_samples_clean %>%
  mutate(analysis.identifier = str_remove(analysis.identifier, pattern_2))

pattern_3 <- or("d", "day") %R% optional("_") %R% capture(one_or_more(DGT))
str_view(P3_all_samples_clean$analysis.identifier, pattern_3)
which_day <- str_match(P3_all_samples_clean$analysis.identifier, pattern_3)

#add day of analysis column
P3_all_samples_clean <- P3_all_samples_clean %>%
  mutate(day = which_day[, 2])

#%>% mutate(analysis.identifier = str_remove(analysis.identifier, pattern_3))

#keep cleaning

pattern_4 <- "COVID"
str_view(P3_all_samples_clean$analysis.identifier, pattern_4)

pattern_5 <- capture("COV" %R% one_or_more(DGT)) %R% "-healthy"
str_view(P3_all_samples_clean$analysis.identifier, pattern_5)

pattern_6 <- or("COV", "COVS") %R% one_or_more(DGT)
str_view(P3_all_samples_clean$analysis.identifier, pattern_6)

pattern_8 <- capture(START %R% WRD %R% WRD %R% optional(DGT)) %R% "_"
pattern_8
str_view(P3_all_samples_clean$analysis.identifier, pattern_8)

pattern_9 <- "_.*"
str_view(P3_all_samples_clean$analysis.identifier, pattern_9)

pattern_10 <- capture("week") %R% "_" %R% capture(DGT)
str_view(P3_all_samples_clean$analysis.identifier, pattern_10)



P3_all_samples_clean <- P3_all_samples_clean %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, pattern_4, "COV")) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, "COV-s", "COVS")) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, "COV-", "COV")) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, pattern_5, REF1)) %>%
  mutate(analysis.identifier = str_replace_all(analysis.identifier, "-", "_")) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, "COVS_", "COVS")) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, "COV0", "COV")) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, " ", "_")) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, pattern_8, str_c(REF1, ".", sep = ""))) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, pattern_10, str_c(REF1, ".", REF2, sep = ""))) %>%
  mutate(analysis.identifier = str_remove(analysis.identifier, pattern_9)) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, "COV12d14", "COV12")) %>%
  mutate(analysis.identifier = str_to_upper(analysis.identifier)) %>%
  mutate(analysis.identifier = replace(analysis.identifier, sample_ID == 	"Specimen_001_ALL_002.fcs", "HD.ALL.week.1")) %>%
  mutate(analysis.identifier = replace(analysis.identifier, sample_ID == 	"Specimen_001_ALL-week-2_002.fcs", "HD.ALL.week.2")) %>%
  mutate(analysis.identifier = replace(analysis.identifier, sample_ID == 	"Specimen_001_ALL-week-3_002.fcs", "HD.ALL.week.3")) %>%
  mutate(analysis.identifier = replace(analysis.identifier, sample_ID == 	"Specimen_001_ALL-week-4_002.fcs", "HD.ALL.week.4")) %>%
  mutate(analysis.identifier = replace(analysis.identifier, analysis.identifier == 	"LF.REPLICATE", "HD.LF.week.1")) %>%
  mutate(analysis.identifier = replace(analysis.identifier, analysis.identifier == 	"SW.REPLICATE", "HD.SW.week.1"))

# for P3 also need to change CK.REPLICATE to HD.CK.week.1

P3_all_samples_clean <- P3_all_samples_clean %>%
  mutate(analysis.identifier = replace(analysis.identifier, analysis.identifier == 	"CK.REPLICATE", "HD.CK.week.1"))
  

pattern_11 <- "HD." %R% capture(one_or_more(WRD))
str_view(P3_all_samples_clean$analysis.identifier, pattern_11)
str_match(P3_all_samples_clean$analysis.identifier, pattern_11)

#count >COV samples with >1 time point and HDs that have been analyzed >1

P3_all_samples_clean %>%
  count(analysis.identifier) %>%
  filter(n > 1)

P3_all_samples_clean %>%
  filter(str_length(analysis.identifier) < 6) 

pattern_12 <- START %R% capture(or("DP", "LW", "MS", "SF", "DB", "GL", "LA", "MT", "LJ", "WA", "NH", "RC", "WM")) %R% END
str_view(P3_all_samples_clean$analysis.identifier, pattern_12)

P3_all_samples_clean$analysis.identifier <- str_replace(P3_all_samples_clean$analysis.identifier, pattern_12, str_c("HD.", REF1, ".week.1"))

pattern_13 <- capture(WRD %R% WRD %R% optional(WRD) %R% ".WEEK." %R% DGT)
str_view(P3_all_samples_clean$analysis.identifier, pattern_13)

P3_all_samples_clean <- P3_all_samples_clean %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, pattern_13, str_c("HD.", REF1, sep = ""))) %>%
  mutate(analysis.identifier = str_to_upper(analysis.identifier)) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, "HD1.T", "HD"))

P3_all_samples_clean <- P3_all_samples_clean %>%
  mutate(analysis.identifier = str_to_upper(analysis.identifier)) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, "HD1.T", "HD"))

pattern_14 <- START %R% "HD.C" %R% END
str_view(P3_all_samples_clean$analysis.identifier, pattern_14)

P3_all_samples_clean <- P3_all_samples_clean %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, pattern_14, "HD.CK")) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, START %R% "HD.R" %R% END, "HD.RC")) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, START %R% "HD.SEPH" %R% END, "HD.SS"))

pattern_15 <- "HD." %R% capture(WRD %R% WRD)
str_view(P3_all_samples_clean$analysis.identifier, pattern_15)


#extract id of HDs

HD_identity <- str_match(P3_all_samples_clean$analysis.identifier, pattern_15)

P3_all_samples_clean <- P3_all_samples_clean %>%
  mutate(HD_id = HD_identity[, 2])

#extract week of time course exp
pattern_16 <- ".WEEK." %R% capture(DGT)
str_view(P3_all_samples_clean$analysis.identifier, pattern_16)


week <- str_match(P3_all_samples_clean$analysis.identifier, pattern_16)

P3_all_samples_clean <- P3_all_samples_clean %>%
  mutate(time_course_week = week[, 2])


#count >COV samples with >1 time point

P3_all_samples_clean %>% filter(subgroup == "COV") %>%
  count(analysis.identifier) %>%
  filter(n > 1)

#add missing time point information

P3_all_samples_clean <- P3_all_samples_clean %>% 
  mutate(day = replace(day, sample_ID == "Specimen_001_COV011 2020-04-22_003.fcs", 9)) %>%
  mutate(day = replace(day, sample_ID == "Specimen_001_COV11-2020-05-06_003.fcs", 24)) %>%
  mutate(day = replace(day, sample_ID == "Specimen_001_COV12-2020-05-06_004.fcs", 22)) %>%
  mutate(day = replace(day, sample_ID == "Specimen_001_COV17_2020-05-02_003.fcs", 3)) %>%
  mutate(day = replace(day, sample_ID == "Specimen_001_COV-18-2020-05-05_005.fcs", 3)) %>%
  mutate(day = replace(day, sample_ID == "Specimen_001_COV-20_2020-06-25_003.fcs", 50)) %>%
  mutate(day = replace(day, sample_ID == "Specimen_001_COV-6_2020-06-02_006.fcs", 61)) %>%
  mutate(day = replace(day, sample_ID == "Specimen_001_COV-s5-2020-05-04_003.fcs", 4)) %>%
  mutate(day = replace(day, sample_ID == "Specimen_001_COV-S7-2020-05-15_004.fcs", 3)) %>%
  filter(analysis.identifier != "IBD686") %>%
  mutate(analysis.identifier = replace(analysis.identifier, sample_ID == "Specimen_001_COV-24_2020-06-02_007.fcs", "COV27"))

pattern_17 <- START %R% capture("COV" %R% optional("S")) %R% capture(DGT) %R% END
str_view(P3_all_samples_clean$analysis.identifier, pattern_17)

pattern_18 <- START %R% capture("COV" %R% optional("S")) %R% capture(DGT %R% DGT) %R% END
str_view(P3_all_samples_clean$analysis.identifier, pattern_18)

P3_all_samples_clean <- P3_all_samples_clean %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, pattern_17, str_c(REF1, "00", REF2))) %>%
  mutate(analysis.identifier = str_replace(analysis.identifier, pattern_18, str_c(REF1, "0", REF2))) %>%
  mutate(day = replace(day, which(is.na(day) & subgroup == "COV"), 0))

#remove technical duplicates


pattern_replicate <- or("Repl-2", "Repl-3", "Replicate-2", "Replicate-3")
str_view(P3_all_samples_clean$sample_ID, pattern_replicate)


P3_all_samples_clean <- P3_all_samples_clean %>%
  filter(!str_detect(sample_ID, pattern_replicate))

#remove template HD-LW
pattern_lw <- "HD-LW"
str_view(P3_all_samples_clean$sample_ID, pattern_lw)

P3_all_samples_clean <- P3_all_samples_clean %>%
  filter(!str_detect(sample_ID, pattern_lw))

rm(list=c('pattern_1','pattern_2','pattern_3','pattern_4','pattern_5','pattern_6','pattern_8',
          'pattern_9','pattern_10','pattern_11','pattern_12','pattern_13','pattern_14','pattern_15',
          'pattern_16','pattern_17','pattern_18','pattern_IBD_SY','pattern_replicate', 'pattern_lw'))

P3_all_samples_clean_join <- P3_all_samples_clean %>% 
  left_join(meta, by = "analysis.identifier") %>% 
  filter(subgroup =="COV")


order_panel_add_ID <- function(x) {
  sg1 <- x[x$subgroup=="COV",]
  sg1_ord <- sg1[order(sg1$analysis.identifier,sg1$day),]
  sg1_ord$sample_order <- NA
  
  names <- levels(as.factor(sg1_ord$analysis.identifier))
  l_names <- length(levels(as.factor(sg1_ord$analysis.identifier)))
  for (i in 1:l_names){
    c <- order(as.integer(sg1_ord$day[sg1_ord$analysis.identifier==names[i]]))
    d <- as.integer(sg1_ord$day[sg1_ord$analysis.identifier==names[i]])
    e <- d[c]
    
    for (j in 1:length(c)){
      sg1_ord$sample_order[sg1_ord$analysis.identifier==names[i] & sg1_ord$day==e[j]] <- j
    }
  }
  sg1_ord$merge_ID <- paste(sg1_ord$analysis.identifier, sg1_ord$sample_order, sep="_")
  
  sg2 <- x[x$subgroup!="COV",]
  sg2$sample_order <- NA
  sg2$merge_ID <- NA
  clean <- rbind(sg1_ord, sg2)
  
  
  return (clean)
}
p3 <- order_panel_add_ID(P3_all_samples_clean_join)



sample_meta <- read_delim("~/Ped-Covid Studie/metadaten/sample_metadata_20210511.csv",
                          ";", escape_double = FALSE, col_types = cols(date_inclusion = col_date(format = "%d.%m.%Y"),
                                                                       date_of_hospital_arrival = col_date(format = "%d.%m.%Y"),
                                                                       sampling_date_ = col_date(format = "%d.%m.%Y")), 
                          trim_ws = TRUE)

sample_meta <- sample_meta[,c(1,2,5:9,26)]
colnames(sample_meta)<- c("sample_ID", "analysis.identifier", "sampling_date", "fever",
                          "WHO_on_date", "medication(0=none,1=NSAID,2=AB,3=IS,4=multi",
                          "time_point", "FACS")

#sample_meta$FACS <- ifelse(sample_meta$analysis.identifier %in% P3_all_samples_clean_join$analysis.identifier, "y", "n")

order_samples_add_ID <- function(x) {
  data <- x[x$FACS=="y",] 
  data$sample_order <- NA
  
  names <- levels(as.factor(data$analysis.identifier))
  
  l_names <- length(levels(as.factor(data$analysis.identifier)))
  
  for (i in 1:l_names){
    c <- order(as.integer(data$time_point[data$analysis.identifier==names[i]]))
    d <- as.integer(data$time_point[data$analysis.identifier==names[i]])
    e <- d[c]
    for (j in 1:length(c)){
      data$sample_order[data$analysis.identifier==names[i] & data$time_point==e[j]] <- j
    }
  }
  
  data$merge_ID <- paste(data$analysis.identifier, data$sample_order, sep="_")
  
  rest <- x[x$FACS != "y",]
  rest$sample_order <- NA
  rest$merge_ID <- NA
  
  data <- rbind(data, rest)
  
  return (data)
}
#order sample meta data
sm <- order_samples_add_ID(sample_meta)

P3_cov_complete <- p3 %>% filter(subgroup=="COV") %>%  left_join(sm, by = "merge_ID")



rm(list=c('HD_identity','p3','P3_all_samples','P3_all_samples_clean','P3_all_samples_clean_join',
          'sm','week','which_day'))
#QC:
P3_cov_complete$analysis.identifier.x[duplicated(P3_cov_complete$merge_ID)]
