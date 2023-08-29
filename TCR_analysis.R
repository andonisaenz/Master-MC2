##################################### IMPLEMENT MIXCR PIPELINE ##################################### 

# This part was executed in the Windows command prompt (Windows Powershell)

# Implement the MiXCR pipeline (rnaseq-cdr3)

# java -jar C:\mixcr\mixcr.jar analyze rnaseq-cdr3 --species mmu 
# C:\Users\andoni\OneDrive\Escritorio\FASTQ\P48\P48_1.fq.gz 
# C:\Users\andoni\OneDrive\Escritorio\FASTQ\P48\P48_2.fq.gz 
# C:\Users\andoni\OneDrive\Escritorio\FASTQ\P48\output_48

# Generate quality control analysis report 

# java -jar C:\mixcr\mixcr.jar exportQc chainUsage 
# C:\Users\andoni\OneDrive\Escritorio\FASTQ\P48\output_48.clns chainUsage_48.pdf

################################# ANALYSIS WITH IMMUNARCH ################################# 

library(immunarch)
library(ggplot2)
library(tidyverse)

path <- "C:\\Users\\andoni\\OneDrive\\Escritorio\\FASTQ"

immdata <- repLoad(path)
immdata_filtered <- immdata$data[!grepl("clones_IG", names(immdata$data))]

CD4_pos <- c(1:3, 7:9, 13:15, 19:21, 25:27, 31:33, 37:39, 43:45, 49:51)
immdata_CD4 <- immdata_filtered[CD4_pos]

CD8_pos <- c(4:6, 10:12, 16:18, 22:24, 28:30, 34:36, 40:42, 46:48, 52:54)
immdata_CD8 <- immdata_filtered[CD8_pos]

#---------------------------------- DATA EXPLORATION ----------------------------------

# Volume

exp_vol_CD4 <- repExplore(immdata_CD4, .method = "volume")
exp_vol_CD8 <- repExplore(immdata_CD8, .method = "volume")

exp_vol_CD4[["Type"]] <- gsub("output_\\d+\\.clones_(TRAD|TRB|TRG)", "\\1", exp_vol_CD4[["Sample"]])
rownames(exp_vol_CD4) <- NULL

exp_vol_CD8[["Type"]] <- gsub("output_\\d+\\.clones_(TRAD|TRB|TRG)", "\\1", exp_vol_CD8[["Sample"]])
rownames(exp_vol_CD8) <- NULL

sample_names_CD4 <- c(rep("P47",3), rep("P49",3), rep("P51",3), rep("P71",3), rep("P73",3), rep("P75",3), 
                      rep("P83",3), rep("P85",3), rep("P87",3))

sample_names_CD8 <- c(rep("P48",3), rep("P50",3), rep("P52",3), rep("P72",3), rep("P74",3), rep("P76",3), 
                      rep("P84",3), rep("P86",3), rep("P88",3))

exp_vol_CD4[["Sample"]] <- sample_names_CD4
exp_vol_CD4 <- exp_vol_CD4[, c("Sample", "Type", "Volume")]
colors <- c(TRAD = "coral", TRB = "dodgerblue", TRG = "springgreen")

exp_vol_CD8[["Sample"]] <- sample_names_CD8
exp_vol_CD8 <- exp_vol_CD8[, c("Sample", "Type", "Volume")]

vol_CD4 <- exp_vol_CD4 %>% 
              ggplot(aes(x = Sample, y = Volume, fill = Type))+ 
              geom_bar(stat = "identity", position = "dodge", width = 0.7)+
              scale_fill_manual(values = colors)+
              labs(x = "Samples", y = "Number of unique clonotypes")+
              theme_minimal()

vol_CD8 <- exp_vol_CD8 %>% 
              ggplot(aes(x = Sample, y = Volume, fill = Type))+ 
              geom_bar(stat = "identity", position = "dodge", width = 0.7)+
              scale_fill_manual(values = colors)+
              labs(x = "Samples", y = "Number of unique clonotypes")+
              theme_minimal()

# Clones

exp_clones_CD4 <- repExplore(immdata_CD4, .method = "clones")
exp_clones_CD8 <- repExplore(immdata_CD8, .method = "clones")

exp_clones_CD4[["Type"]] <- gsub("output_\\d+\\.clones_(TRAD|TRB|TRG)", "\\1", exp_clones_CD4[["Sample"]])
rownames(exp_clones_CD4) <- NULL
exp_clones_CD4[["Sample"]] <- sample_names_CD4
exp_clones_CD4 <- exp_clones_CD4[, c("Sample", "Type", "Clones")]

exp_clones_CD8[["Type"]] <- gsub("output_\\d+\\.clones_(TRAD|TRB|TRG)", "\\1", exp_clones_CD8[["Sample"]])
rownames(exp_clones_CD8) <- NULL
exp_clones_CD8[["Sample"]] <- sample_names_CD8
exp_clones_CD8 <- exp_clones_CD8[, c("Sample", "Type", "Clones")]


clones_CD4 <- exp_clones_CD4 %>% 
                ggplot(aes(x = Sample, y = Clones, fill = Type, width = 0.8))+
                geom_col()+
                scale_fill_manual(values = colors)+
                labs(x = "Samples", y = "Number of total reads")+
                theme_minimal()

clones_CD8 <- exp_clones_CD8 %>% 
                ggplot(aes(x = Sample, y = Clones, fill = Type, width = 0.8))+
                geom_col()+
                scale_fill_manual(values = colors)+
                labs(x = "Samples", y = "Number of total reads")+
                theme_minimal()

ggarrange(vol_CD4, vol_CD8, clones_CD4, clones_CD8, ncol = 2, nrow = 2, common.legend = TRUE, 
          legend = "bottom", labels = c("A", "B", "C", "D"))

# Length

exp_len_CD4 <- repExplore(immdata_CD4, .method = "len", .col = "aa")
exp_len_CD4[["Sample"]] <- gsub("output_(\\d+)\\.clones_(\\w+)", "P\\1_\\2", exp_len_CD4[["Sample"]])

exp_len_CD4_untreated <- exp_len_CD4[grep("P47", exp_len_CD4$Sample), ]
exp_len_CD4_treated <- exp_len_CD4[grep("P71", exp_len_CD4$Sample), ]
exp_len_CD4_control <- exp_len_CD4[grep("P83", exp_len_CD4$Sample), ]

len_CD4_untreated <- vis(exp_len_CD4_untreated)
len_CD4_treated <- vis(exp_len_CD4_treated)
len_CD4_control <- vis(exp_len_CD4_control)

exp_len_CD8 <- repExplore(immdata_CD8, .method = "len", .col = "aa")
exp_len_CD8[["Sample"]] <- gsub("output_(\\d+)\\.clones_(\\w+)", "P\\1_\\2", exp_len_CD8[["Sample"]])

exp_len_CD8_untreated <- exp_len_CD8[grep("P48", exp_len_CD8$Sample), ]
exp_len_CD8_treated <- exp_len_CD8[grep("P72", exp_len_CD8$Sample), ]
exp_len_CD8_control <- exp_len_CD8[grep("P84", exp_len_CD8$Sample), ]

len_CD8_untreated <- vis(exp_len_CD8_untreated)
len_CD8_treated <- vis(exp_len_CD8_treated)
len_CD8_control <- vis(exp_len_CD8_control)

ggarrange(len_CD4_untreated, len_CD4_treated, len_CD4_control, 
          len_CD8_untreated, len_CD8_treated, len_CD8_control,
          ncol = 3, nrow = 2, labels = c("A", "B", "C", "D", "E", "F"))

# Counts

exp_count_CD4 <- repExplore(immdata_CD4, .method = "count")
exp_count_CD4[["Sample"]] <- gsub("output_(\\d+)\\.clones_(\\w+)", "P\\1_\\2", exp_count_CD4[["Sample"]])

exp_count_CD4_untreated <- exp_count_CD4[grep("P47", exp_count_CD4$Sample), ]
exp_count_CD4_treated <- exp_count_CD4[grep("P71", exp_count_CD4$Sample), ]
exp_count_CD4_control <- exp_count_CD4[grep("P83", exp_count_CD4$Sample), ]

count_CD4_untreated <- vis(exp_count_CD4_untreated)
count_CD4_treated <- vis(exp_count_CD4_treated)
count_CD4_control <- vis(exp_count_CD4_control)

exp_count_CD8 <- repExplore(immdata_CD8, .method = "count")
exp_count_CD8[["Sample"]] <- gsub("output_(\\d+)\\.clones_(\\w+)", "P\\1_\\2", exp_count_CD8[["Sample"]])

exp_count_CD8_untreated <- exp_count_CD8[grep("P48", exp_count_CD8$Sample), ]
exp_count_CD8_treated <- exp_count_CD8[grep("P72", exp_count_CD8$Sample), ]
exp_count_CD8_control <- exp_count_CD8[grep("P84", exp_count_CD8$Sample), ]

count_CD8_untreated <- vis(exp_count_CD8_untreated)
count_CD8_treated <- vis(exp_count_CD8_treated)
count_CD8_control <- vis(exp_count_CD8_control)

ggarrange(count_CD4_untreated, count_CD4_treated, count_CD4_control, 
          count_CD8_untreated, count_CD8_treated, count_CD8_control,
          ncol = 3, nrow = 2, labels = c("A", "B", "C", "D", "E", "F"))

#---------------------------------- REPERTOIRE ANALYSIS ----------------------------------

library(packcircles)
library(ggforce)
library(ggpubr)

# CD4 untreated

immdata_CD4_untreated_TRAD <- rbind(immdata_CD4$output_47.clones_TRAD, 
                             immdata_CD4$output_49.clones_TRAD, 
                             immdata_CD4$output_51.clones_TRAD)
immdata_CD4_untreated_TRAD <- immdata_CD4_untreated_TRAD[order(immdata_CD4_untreated_TRAD$Proportion, decreasing = T),]

immdata_CD4_untreated_TRB <- rbind(immdata_CD4$output_47.clones_TRB, 
                                    immdata_CD4$output_49.clones_TRB, 
                                    immdata_CD4$output_51.clones_TRB)
immdata_CD4_untreated_TRB <- immdata_CD4_untreated_TRB[order(immdata_CD4_untreated_TRB$Proportion, decreasing = T),]

immdata_CD4_untreated_TRG <- rbind(immdata_CD4$output_47.clones_TRG, 
                                   immdata_CD4$output_49.clones_TRG, 
                                   immdata_CD4$output_51.clones_TRG)
immdata_CD4_untreated_TRG <- immdata_CD4_untreated_TRG[order(immdata_CD4_untreated_TRG$Proportion, decreasing = T),]

total_clones_CD4_untreated <- sum(immdata_CD4_untreated_TRAD$Clones, 
                                  immdata_CD4_untreated_TRB$Clones, 
                                  immdata_CD4_untreated_TRG$Clones)

immdata_CD4_untreated_TRAD$Prop <- 100*(immdata_CD4_untreated_TRAD$Clones/total_clones_CD4_untreated)
immdata_CD4_untreated_TRB$Prop <- 100*(immdata_CD4_untreated_TRB$Clones/total_clones_CD4_untreated)
immdata_CD4_untreated_TRG$Prop <- 100*(immdata_CD4_untreated_TRG$Clones/total_clones_CD4_untreated)

packing_CD4_untreated_TRAD <- circleProgressiveLayout(immdata_CD4_untreated_TRAD$Prop, sizetype = "area")
packing_CD4_untreated_TRB <- circleProgressiveLayout(immdata_CD4_untreated_TRB$Prop, sizetype = "area")
packing_CD4_untreated_TRG <- circleProgressiveLayout(immdata_CD4_untreated_TRG$Prop, sizetype = "area")

packing_CD4_untreated_TRAD <- cbind(immdata_CD4_untreated_TRAD, packing_CD4_untreated_TRAD)
packing_CD4_untreated_TRB <- cbind(immdata_CD4_untreated_TRB, packing_CD4_untreated_TRB)
packing_CD4_untreated_TRG <- cbind(immdata_CD4_untreated_TRG, packing_CD4_untreated_TRG)

packing_CD4_untreated_TRAD$x <- packing_CD4_untreated_TRAD$x - 10
packing_CD4_untreated_TRB$x <- packing_CD4_untreated_TRB$x + 10
packing_CD4_untreated_TRG$y <- packing_CD4_untreated_TRG$y - 10

packing_CD4_untreated <- rbind(packing_CD4_untreated_TRAD, packing_CD4_untreated_TRB, packing_CD4_untreated_TRG)

f1 <- ggplot() + 
  geom_circle(data = packing_CD4_untreated, aes(x0 = x, y0 = y, r = radius, fill = cut(Clones, c(0,10,50,100,300, max(Clones)))), linetype = 0)+
  scale_fill_manual(values = c("grey", "yellow", "orange", "red", "darkorchid"), 
                    labels = c("< 10", "10-50", "50-100", "100-300", "> 300"))+
  geom_text(data = subset(packing_CD4_untreated, Clones > 100), aes(x, y, size = Clones, 
                                                                   label = paste(str_extract(V.name, "(TRA|TRB|TRG)V\\d{1,2}"), 
                                                                                 "\n", 
                                                                                 str_extract(J.name, "(TRA|TRB|TRG)J\\d{1,2}"))))+
  scale_size_continuous(range = c(1.5,2.5))+
  theme_void()+ 
  theme(legend.position = "right")+
  coord_equal()+
  labs(title = "Untreated · CD4")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(size = "none")+
  labs(fill = "Clone size")

# CD4 treated

immdata_CD4_treated_TRAD <- rbind(immdata_CD4$output_71.clones_TRAD, 
                                    immdata_CD4$output_73.clones_TRAD, 
                                    immdata_CD4$output_75.clones_TRAD)
immdata_CD4_treated_TRAD <- immdata_CD4_treated_TRAD[order(immdata_CD4_treated_TRAD$Clones, decreasing = T),]

immdata_CD4_treated_TRB <- rbind(immdata_CD4$output_71.clones_TRB, 
                                   immdata_CD4$output_73.clones_TRB, 
                                   immdata_CD4$output_75.clones_TRB)
immdata_CD4_treated_TRB <- immdata_CD4_treated_TRB[order(immdata_CD4_treated_TRB$Clones, decreasing = T),]

immdata_CD4_treated_TRG <- rbind(immdata_CD4$output_71.clones_TRG, 
                                   immdata_CD4$output_73.clones_TRG, 
                                   immdata_CD4$output_75.clones_TRG)
immdata_CD4_treated_TRG <- immdata_CD4_treated_TRG[order(immdata_CD4_treated_TRG$Clones, decreasing = T),]

total_clones_CD4_treated <- sum(immdata_CD4_treated_TRAD$Clones, 
                                immdata_CD4_treated_TRB$Clones, 
                                immdata_CD4_treated_TRG$Clones)

immdata_CD4_treated_TRAD$Prop <- 100*(immdata_CD4_treated_TRAD$Clones/total_clones_CD4_treated)
immdata_CD4_treated_TRB$Prop <- 100*(immdata_CD4_treated_TRB$Clones/total_clones_CD4_treated)
immdata_CD4_treated_TRG$Prop <- 100*(immdata_CD4_treated_TRG$Clones/total_clones_CD4_treated)

packing_CD4_treated_TRAD <- circleProgressiveLayout(immdata_CD4_treated_TRAD$Prop, sizetype = "area")
packing_CD4_treated_TRB <- circleProgressiveLayout(immdata_CD4_treated_TRB$Prop, sizetype = "area")
packing_CD4_treated_TRG <- circleProgressiveLayout(immdata_CD4_treated_TRG$Prop, sizetype = "area")

packing_CD4_treated_TRAD <- cbind(immdata_CD4_treated_TRAD, packing_CD4_treated_TRAD)
packing_CD4_treated_TRB <- cbind(immdata_CD4_treated_TRB, packing_CD4_treated_TRB)
packing_CD4_treated_TRG <- cbind(immdata_CD4_treated_TRG, packing_CD4_treated_TRG)

packing_CD4_treated_TRAD$x <- packing_CD4_treated_TRAD$x - 10
packing_CD4_treated_TRB$x <- packing_CD4_treated_TRB$x + 10
packing_CD4_treated_TRG$y <- packing_CD4_treated_TRG$y - 10

packing_CD4_treated <- rbind(packing_CD4_treated_TRAD, packing_CD4_treated_TRB, packing_CD4_treated_TRG)

f2 <- ggplot() + 
  geom_circle(data = packing_CD4_treated, aes(x0 = x, y0 = y, r = radius, fill = cut(Clones, c(0,10,50,100,300, max(Clones)))), linetype = 0)+
  scale_fill_manual(values = c("grey", "yellow", "orange", "red", "darkorchid"), 
                    labels = c("< 10", "10-50", "50-100", "100-300", "> 300"))+
  geom_text(data = subset(packing_CD4_treated, Clones > 100), aes(x, y, size = Clones, 
                                                                    label = paste(str_extract(V.name, "(TRA|TRB|TRG)V\\d{1,2}"), 
                                                                                  "\n", 
                                                                                  str_extract(J.name, "(TRA|TRB|TRG)J\\d{1,2}"))))+
  scale_size_continuous(range = c(1.5,2.5))+
  theme_void()+ 
  theme(legend.position = "right")+
  coord_equal()+
  labs(title = "Treated · CD4")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(size = "none")+
  labs(fill = "Clone size")

# CD4 control

immdata_CD4_control_TRAD <- rbind(immdata_CD4$output_83.clones_TRAD, 
                                  immdata_CD4$output_85.clones_TRAD, 
                                  immdata_CD4$output_87.clones_TRAD)
immdata_CD4_control_TRAD <- immdata_CD4_control_TRAD[order(immdata_CD4_control_TRAD$Clones, decreasing = T),]

immdata_CD4_control_TRB <- rbind(immdata_CD4$output_83.clones_TRB, 
                                 immdata_CD4$output_85.clones_TRB, 
                                 immdata_CD4$output_87.clones_TRB)
immdata_CD4_control_TRB <- immdata_CD4_control_TRB[order(immdata_CD4_control_TRB$Clones, decreasing = T),]

immdata_CD4_control_TRG <- rbind(immdata_CD4$output_83.clones_TRG, 
                                 immdata_CD4$output_85.clones_TRG, 
                                 immdata_CD4$output_87.clones_TRG)
immdata_CD4_control_TRG <- immdata_CD4_control_TRG[order(immdata_CD4_control_TRG$Clones, decreasing = T),]

total_clones_CD4_control <- sum(immdata_CD4_control_TRAD$Clones, 
                                immdata_CD4_control_TRB$Clones, 
                                immdata_CD4_control_TRG$Clones)

immdata_CD4_control_TRAD$Prop <- 100*(immdata_CD4_control_TRAD$Clones/total_clones_CD4_control)
immdata_CD4_control_TRB$Prop <- 100*(immdata_CD4_control_TRB$Clones/total_clones_CD4_control)
immdata_CD4_control_TRG$Prop <- 100*(immdata_CD4_control_TRG$Clones/total_clones_CD4_control)

packing_CD4_control_TRAD <- circleProgressiveLayout(immdata_CD4_control_TRAD$Prop, sizetype = "area")
packing_CD4_control_TRB <- circleProgressiveLayout(immdata_CD4_control_TRB$Prop, sizetype = "area")
packing_CD4_control_TRG <- circleProgressiveLayout(immdata_CD4_control_TRG$Prop, sizetype = "area")

packing_CD4_control_TRAD <- cbind(immdata_CD4_control_TRAD, packing_CD4_control_TRAD)
packing_CD4_control_TRB <- cbind(immdata_CD4_control_TRB, packing_CD4_control_TRB)
packing_CD4_control_TRG <- cbind(immdata_CD4_control_TRG, packing_CD4_control_TRG)

packing_CD4_control_TRAD$x <- packing_CD4_control_TRAD$x - 10
packing_CD4_control_TRB$x <- packing_CD4_control_TRB$x + 10
packing_CD4_control_TRG$y <- packing_CD4_control_TRG$y - 10

packing_CD4_control <- rbind(packing_CD4_control_TRAD, packing_CD4_control_TRB, packing_CD4_control_TRG)

f3 <- ggplot() + 
  geom_circle(data = packing_CD4_control, aes(x0 = x, y0 = y, r = radius, fill = cut(Clones, c(0,10,50,100,300, max(Clones)))), linetype = 0)+
  scale_fill_manual(values = c("grey", "yellow", "orange", "red", "darkorchid"), 
                    labels = c("< 10", "10-50", "50-100", "100-300", "> 300"))+
  geom_text(data = subset(packing_CD4_control, Clones > 100), aes(x, y, size = Clones, 
                                                                  label = paste(str_extract(V.name, "(TRA|TRB|TRG)V\\d{1,2}"), 
                                                                                "\n", 
                                                                                str_extract(J.name, "(TRA|TRB|TRG)J\\d{1,2}"))))+
  scale_size_continuous(range = c(1.5,2.5))+
  theme_void()+ 
  theme(legend.position = "right")+
  coord_equal()+
  labs(title = "Control · CD4")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(size = "none")+
  labs(fill = "Clone size")

# CD8 untreated

immdata_CD8_untreated_TRAD <- rbind(immdata_CD8$output_48.clones_TRAD, 
                                    immdata_CD8$output_50.clones_TRAD, 
                                    immdata_CD8$output_52.clones_TRAD)
immdata_CD8_untreated_TRAD <- immdata_CD8_untreated_TRAD[order(immdata_CD8_untreated_TRAD$Clones, decreasing = T),]

immdata_CD8_untreated_TRB <- rbind(immdata_CD8$output_48.clones_TRB, 
                                   immdata_CD8$output_50.clones_TRB, 
                                   immdata_CD8$output_52.clones_TRB)
immdata_CD8_untreated_TRB <- immdata_CD8_untreated_TRB[order(immdata_CD8_untreated_TRB$Clones, decreasing = T),]

immdata_CD8_untreated_TRG <- rbind(immdata_CD8$output_48.clones_TRG, 
                                   immdata_CD8$output_50.clones_TRG, 
                                   immdata_CD8$output_52.clones_TRG)
immdata_CD8_untreated_TRG <- immdata_CD8_untreated_TRG[order(immdata_CD8_untreated_TRG$Clones, decreasing = T),]

total_clones_CD8_untreated <- sum(immdata_CD8_untreated_TRAD$Clones, 
                                  immdata_CD8_untreated_TRB$Clones, 
                                  immdata_CD8_untreated_TRG$Clones)

immdata_CD8_untreated_TRAD$Prop <- 100*(immdata_CD8_untreated_TRAD$Clones/total_clones_CD8_untreated)
immdata_CD8_untreated_TRB$Prop <- 100*(immdata_CD8_untreated_TRB$Clones/total_clones_CD8_untreated)
immdata_CD8_untreated_TRG$Prop <- 100*(immdata_CD8_untreated_TRG$Clones/total_clones_CD8_untreated)

packing_CD8_untreated_TRAD <- circleProgressiveLayout(immdata_CD8_untreated_TRAD$Prop, sizetype = "area")
packing_CD8_untreated_TRB <- circleProgressiveLayout(immdata_CD8_untreated_TRB$Prop, sizetype = "area")
packing_CD8_untreated_TRG <- circleProgressiveLayout(immdata_CD8_untreated_TRG$Prop, sizetype = "area")

packing_CD8_untreated_TRAD <- cbind(immdata_CD8_untreated_TRAD, packing_CD8_untreated_TRAD)
packing_CD8_untreated_TRB <- cbind(immdata_CD8_untreated_TRB, packing_CD8_untreated_TRB)
packing_CD8_untreated_TRG <- cbind(immdata_CD8_untreated_TRG, packing_CD8_untreated_TRG)

packing_CD8_untreated_TRAD$x <- packing_CD8_untreated_TRAD$x - 10
packing_CD8_untreated_TRB$x <- packing_CD8_untreated_TRB$x + 10
packing_CD8_untreated_TRG$y <- packing_CD8_untreated_TRG$y - 10

packing_CD8_untreated <- rbind(packing_CD8_untreated_TRAD, packing_CD8_untreated_TRB, packing_CD8_untreated_TRG)

f4 <- ggplot() + 
  geom_circle(data = packing_CD8_untreated, aes(x0 = x, y0 = y, r = radius, fill = cut(Clones, c(0,10,50,100,300, max(Clones)))), linetype = 0)+
  scale_fill_manual(values = c("grey", "yellow", "orange", "red", "darkorchid"), 
                    labels = c("< 10", "10-50", "50-100", "100-300", "> 300"))+
  geom_text(data = subset(packing_CD8_untreated, Clones > 100), aes(x, y, size = Clones, 
                                                                  label = paste(str_extract(V.name, "(TRA|TRB|TRG)V\\d{1,2}"), 
                                                                                "\n", 
                                                                                str_extract(J.name, "(TRA|TRB|TRG)J\\d{1,2}"))))+
  scale_size_continuous(range = c(1.5,2.5))+
  theme_void()+ 
  theme(legend.position = "right")+
  coord_equal()+
  labs(title = "Untreated · CD8")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(size = "none")+
  labs(fill = "Clone size")

# CD8 treated

immdata_CD8_treated_TRAD <- rbind(immdata_CD8$output_72.clones_TRAD, 
                                  immdata_CD8$output_74.clones_TRAD, 
                                  immdata_CD8$output_76.clones_TRAD)
immdata_CD8_treated_TRAD <- immdata_CD8_treated_TRAD[order(immdata_CD8_treated_TRAD$Clones, decreasing = T),]

immdata_CD8_treated_TRB <- rbind(immdata_CD8$output_72.clones_TRB, 
                                 immdata_CD8$output_74.clones_TRB, 
                                 immdata_CD8$output_76.clones_TRB)
immdata_CD8_treated_TRB <- immdata_CD8_treated_TRB[order(immdata_CD8_treated_TRB$Clones, decreasing = T),]

immdata_CD8_treated_TRG <- rbind(immdata_CD8$output_72.clones_TRG, 
                                 immdata_CD8$output_74.clones_TRG, 
                                 immdata_CD8$output_76.clones_TRG)
immdata_CD8_treated_TRG <- immdata_CD8_treated_TRG[order(immdata_CD8_treated_TRG$Clones, decreasing = T),]

total_clones_CD8_treated <- sum(immdata_CD8_treated_TRAD$Clones, 
                                immdata_CD8_treated_TRB$Clones, 
                                immdata_CD8_treated_TRG$Clones)

immdata_CD8_treated_TRAD$Prop <- 100*(immdata_CD8_treated_TRAD$Clones/total_clones_CD8_treated)
immdata_CD8_treated_TRB$Prop <- 100*(immdata_CD8_treated_TRB$Clones/total_clones_CD8_treated)
immdata_CD8_treated_TRG$Prop <- 100*(immdata_CD8_treated_TRG$Clones/total_clones_CD8_treated)

packing_CD8_treated_TRAD <- circleProgressiveLayout(immdata_CD8_treated_TRAD$Prop, sizetype = "area")
packing_CD8_treated_TRB <- circleProgressiveLayout(immdata_CD8_treated_TRB$Prop, sizetype = "area")
packing_CD8_treated_TRG <- circleProgressiveLayout(immdata_CD8_treated_TRG$Prop, sizetype = "area")

packing_CD8_treated_TRAD <- cbind(immdata_CD8_treated_TRAD, packing_CD8_treated_TRAD)
packing_CD8_treated_TRB <- cbind(immdata_CD8_treated_TRB, packing_CD8_treated_TRB)
packing_CD8_treated_TRG <- cbind(immdata_CD8_treated_TRG, packing_CD8_treated_TRG)

packing_CD8_treated_TRAD$x <- packing_CD8_treated_TRAD$x - 10
packing_CD8_treated_TRB$x <- packing_CD8_treated_TRB$x + 10
packing_CD8_treated_TRG$y <- packing_CD8_treated_TRG$y - 10

packing_CD8_treated <- rbind(packing_CD8_treated_TRAD, packing_CD8_treated_TRB, packing_CD8_treated_TRG)

f5 <- ggplot() + 
  geom_circle(data = packing_CD8_treated, aes(x0 = x, y0 = y, r = radius, fill = cut(Clones, c(0,10,50,100,300, max(Clones)))), linetype = 0)+
  scale_fill_manual(values = c("grey", "yellow", "orange", "red", "darkorchid"), 
                    labels = c("< 10", "10-50", "50-100", "100-300", "> 300"))+
  geom_text(data = subset(packing_CD8_treated, Clones > 100), aes(x, y, size = Clones, 
                                                                    label = paste(str_extract(V.name, "(TRA|TRB|TRG)V\\d{1,2}"), 
                                                                                  "\n", 
                                                                                  str_extract(J.name, "(TRA|TRB|TRG)J\\d{1,2}"))))+
  scale_size_continuous(range = c(1.5,2.5))+
  theme_void()+ 
  theme(legend.position = "right")+
  coord_equal()+
  labs(title = "Treated · CD8")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(size = "none")+
  labs(fill = "Clone size")

# CD8 control

immdata_CD8_control_TRAD <- rbind(immdata_CD8$output_84.clones_TRAD, 
                                  immdata_CD8$output_86.clones_TRAD, 
                                  immdata_CD8$output_88.clones_TRAD)
immdata_CD8_control_TRAD <- immdata_CD8_control_TRAD[order(immdata_CD8_control_TRAD$Clones, decreasing = T),]

immdata_CD8_control_TRB <- rbind(immdata_CD8$output_84.clones_TRB, 
                                 immdata_CD8$output_86.clones_TRB, 
                                 immdata_CD8$output_88.clones_TRB)
immdata_CD8_control_TRB <- immdata_CD8_control_TRB[order(immdata_CD8_control_TRB$Clones, decreasing = T),]

immdata_CD8_control_TRG <- rbind(immdata_CD8$output_84.clones_TRG, 
                                 immdata_CD8$output_86.clones_TRG, 
                                 immdata_CD8$output_88.clones_TRG)
immdata_CD8_control_TRG <- immdata_CD8_control_TRG[order(immdata_CD8_control_TRG$Clones, decreasing = T),]

total_clones_CD8_control <- sum(immdata_CD8_control_TRAD$Clones, 
                                immdata_CD8_control_TRB$Clones, 
                                immdata_CD8_control_TRG$Clones)

immdata_CD8_control_TRAD$Prop <- 100*(immdata_CD8_control_TRAD$Clones/total_clones_CD8_control)
immdata_CD8_control_TRB$Prop <- 100*(immdata_CD8_control_TRB$Clones/total_clones_CD8_control)
immdata_CD8_control_TRG$Prop <- 100*(immdata_CD8_control_TRG$Clones/total_clones_CD8_control)

packing_CD8_control_TRAD <- circleProgressiveLayout(immdata_CD8_control_TRAD$Prop, sizetype = "area")
packing_CD8_control_TRB <- circleProgressiveLayout(immdata_CD8_control_TRB$Prop, sizetype = "area")
packing_CD8_control_TRG <- circleProgressiveLayout(immdata_CD8_control_TRG$Prop, sizetype = "area")

packing_CD8_control_TRAD <- cbind(immdata_CD8_control_TRAD, packing_CD8_control_TRAD)
packing_CD8_control_TRB <- cbind(immdata_CD8_control_TRB, packing_CD8_control_TRB)
packing_CD8_control_TRG <- cbind(immdata_CD8_control_TRG, packing_CD8_control_TRG)

packing_CD8_control_TRAD$x <- packing_CD8_control_TRAD$x - 10
packing_CD8_control_TRB$x <- packing_CD8_control_TRB$x + 10
packing_CD8_control_TRG$y <- packing_CD8_control_TRG$y - 10

packing_CD8_control <- rbind(packing_CD8_control_TRAD, packing_CD8_control_TRB, packing_CD8_control_TRG)

f6 <- ggplot() + 
  geom_circle(data = packing_CD8_control, aes(x0 = x, y0 = y, r = radius, fill = cut(Clones, c(0,10,50,100,300, max(Clones)))), linetype = 0)+
  scale_fill_manual(values = c("grey", "yellow", "orange", "red", "darkorchid"), 
                    labels = c("< 10", "10-50", "50-100", "100-300", "> 300"))+
  geom_text(data = subset(packing_CD8_control, Clones > 100), aes(x, y, size = Clones, 
                                                                  label = paste(str_extract(V.name, "(TRA|TRB|TRG)V\\d{1,2}"), 
                                                                                "\n", 
                                                                                str_extract(J.name, "(TRA|TRB|TRG)J\\d{1,2}"))))+
  scale_size_continuous(range = c(1.5,2.5))+
  theme_void()+ 
  theme(legend.position = "right")+
  coord_equal()+
  labs(title = "Control · CD8")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(size = "none")+
  labs(fill = "Clone size")

ggarrange(f1, f4, f2, f5, f3, f6, ncol = 2, nrow = 3, common.legend = T,
          legend = "bottom", labels = c("A", "B", "C", "D", "E", "F"))

# ---------------------------------- SHANNON-WIENER INDEX ----------------------------------

library(vegan)
library(ggpubr)

# Get the number of unique clonotypes in each sample type

diversity_clones_CD4 <- data.frame(matrix(data = NA, nrow = 27, ncol = 2, byrow = T))
colnames(diversity_clones_CD4) <- c("Sample", "Type")
diversity_clones_CD4["Sample"] <- sample_names_CD4
diversity_clones_CD4["Type"] <- rep(c("TRAD", "TRB", "TRG"))
treatment <- c(rep(c("Untreated"), 9), rep(c("Treated"), 9), rep(c("Control"), 9))
diversity_clones_CD4["Treatment"] <- treatment
diversity_clones_CD4["Shannon_index"] <- NA

for(i in 1:length(immdata_CD4)){
  diversity_clones_CD4[i,4] <- diversity(immdata_CD4[[i]]$Clones)
}

diversity_clones_CD8 <- data.frame(matrix(data = NA, nrow = 27, ncol = 2, byrow = T))
colnames(diversity_clones_CD8) <- c("Sample", "Type")
diversity_clones_CD8["Sample"] <- sample_names_CD8
diversity_clones_CD8["Type"] <- rep(c("TRAD", "TRB", "TRG"))
diversity_clones_CD8["Treatment"] <- treatment
diversity_clones_CD8["Shannon_index"] <- NA

for(i in 1:length(immdata_CD8)){
  diversity_clones_CD8[i,4] <- diversity(immdata_CD8[[i]]$Clones)
}

sh1 <- diversity_clones_CD4%>% 
  ggplot(aes(x = Treatment, y = Shannon_index, fill = Type))+
  geom_boxplot(alpha = 0.2)+
  labs(y = "Shannon index")+
  ylim(0,8)

# Perform ANOVA tests for TRAD, TRB and TRG (significant difference between group means) in CD4 samples

compare_means(Shannon_index ~ Treatment, data = diversity_clones_CD4[diversity_clones_CD4$Type == "TRAD",], method = "anova")
compare_means(Shannon_index ~ Treatment, data = diversity_clones_CD4[diversity_clones_CD4$Type == "TRB",], method = "anova")
compare_means(Shannon_index ~ Treatment, data = diversity_clones_CD4[diversity_clones_CD4$Type == "TRG",], method = "anova")

# Perform post-hoc t-tests (pairwise comparisons between 2 experimental groups) 

compare_means(Shannon_index ~ Treatment, data = diversity_clones_CD4[diversity_clones_CD4$Type == "TRAD",], method = "t.test")
compare_means(Shannon_index ~ Treatment, data = diversity_clones_CD4[diversity_clones_CD4$Type == "TRB",], method = "t.test")
compare_means(Shannon_index ~ Treatment, data = diversity_clones_CD4[diversity_clones_CD4$Type == "TRG",], method = "t.test")

sh2 <- diversity_clones_CD8 %>% 
  ggplot(aes(x = Treatment, y = Shannon_index, fill = Type))+
  geom_boxplot(alpha = 0.2)+
  labs(y = "Shannon index")+
  ylim(0,8)

# Perform ANOVA tests for TRAD, TRB and TRG in CD8 samples

compare_means(Shannon_index ~ Treatment, data = diversity_clones_CD8[diversity_clones_CD8$Type == "TRAD",], method = "anova")
compare_means(Shannon_index ~ Treatment, data = diversity_clones_CD8[diversity_clones_CD8$Type == "TRB",], method = "anova")
compare_means(Shannon_index ~ Treatment, data = diversity_clones_CD8[diversity_clones_CD8$Type == "TRG",], method = "anova")

# Perform post-hoc t-tests  

compare_means(Shannon_index ~ Treatment, data = diversity_clones_CD8[diversity_clones_CD8$Type == "TRAD",], method = "t.test")
compare_means(Shannon_index ~ Treatment, data = diversity_clones_CD8[diversity_clones_CD8$Type == "TRB",], method = "t.test")
compare_means(Shannon_index ~ Treatment, data = diversity_clones_CD8[diversity_clones_CD8$Type == "TRG",], method = "t.test")

ggarrange(sh1, sh2, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom", labels = c("A", "B"))




