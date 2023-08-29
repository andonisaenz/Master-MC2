
library(readr)
library(tidyr)
library(tidyverse)
library(ggplot2)

MultiQC_statistics <- read_delim("C:/Users/andon/OneDrive/Escritorio/MultiQC_statistics.csv", 
                                 delim = ";", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

sample_names <- as.matrix(MultiQC_statistics[,1])
MultiQC_statistics <- MultiQC_statistics[,-1]
rownames(MultiQC_statistics) <- sample_names

datos <- as.factor(as.matrix(MultiQC_statistics))
samples <- rownames(MultiQC_statistics)
checkpoints <- colnames(MultiQC_statistics)
colores <- c("PASS" = "lightgreen", "WARNING" = "yellow", "ERROR" = "red")

df <- data.frame(Samples = rep(samples, length(checkpoints)),
                 Checkpoints = rep(checkpoints, each = length(samples)),
                 Values = as.vector(datos))

ggplot(df, aes(x = Samples, y = Checkpoints, fill = Values))+
        geom_tile(color = "white")+
        scale_fill_manual(values = colores)+
        theme_minimal()+
        labs(x = "Samples", y = "Checkpoints")+
        coord_fixed(ratio = 1)+
        theme(axis.text.x = element_text(size = 9, angle = 40))
  

#########################################################################################

library(readxl)
library(ggpubr)

Per_base_quality_CD4 <- read_excel("C:/Users/andon/OneDrive/Escritorio/Per_base_sequence_quality_CD4.xlsx")

Per_base_quality_CD8 <- read_excel("C:/Users/andon/OneDrive/Escritorio/Per_base_sequence_quality_CD8.xlsx")

seq_duplication_CD4 <- read_excel("C:/Users/andon/OneDrive/Escritorio/Sequence_duplication_levels_CD4.xlsx")
seq_duplication_CD4 <- pivot_longer(seq_duplication_CD4, cols = 2:19, names_to = "Sample")
seq_duplication_CD4$Level <- as.factor(seq_duplication_CD4$Level)
seq_duplication_CD4$Level <- ordered(seq_duplication_CD4$Level, levels = c("1", "2", "3", "4", "5", "6", "7",
                                                                          "8", "9", ">10", ">50", ">100", 
                                                                          ">500", ">1k", ">5k", ">10k+"))

seq_duplication_CD8 <- read_excel("C:/Users/andon/OneDrive/Escritorio/Sequence_duplication_levels_CD8.xlsx")
seq_duplication_CD8 <- pivot_longer(seq_duplication_CD8, cols = 2:19, names_to = "Sample")
seq_duplication_CD8$Level <- as.factor(seq_duplication_CD8$Level)
seq_duplication_CD8$Level <- ordered(seq_duplication_CD8$Level, levels = c("1", "2", "3", "4", "5", "6", "7",
                                                                           "8", "9", ">10", ">50", ">100", 
                                                                           ">500", ">1k", ">5k", ">10k+"))

p1 <- per_base_GC_CD4 %>% 
  ggplot(aes(x = Base, y = value, color = Sample))+
  geom_line()+
  labs(x = "% GC", y = "Number of reads")

p2 <- per_base_GC_CD8 %>% 
  ggplot(aes(x = Base, y = value, color = Sample))+
  geom_line()+
  labs(x = "% GC", y = "Number of reads")

p3 <- per_base_quality_CD4 %>% 
  ggplot(aes(x = Base, y = value, color = Sample))+
  geom_line()+
  labs(x = "Position (bp)", y = "Phred score")

p4 <- per_base_quality_CD8 %>% 
  ggplot(aes(x = Base, y = value, color = Sample))+
  geom_line()+
  labs(x = "Position (bp)", y = "Phred score")

p5 <- seq_duplication_CD4 %>% 
  ggplot(aes(x = as.numeric(Level), y = value, color = Sample))+
  geom_line()+
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7",
                              "8", "9", ">10", ">50", ">100", 
                              ">500", ">1k", ">5k", ">10k+"))+
  labs(x = "Sequence duplication levels", y = "% of library")

p6 <-seq_duplication_CD8 %>% 
  ggplot(aes(x = as.numeric(Level), y = value, color = Sample))+
  geom_line()+
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7",
                              "8", "9", ">10", ">50", ">100", 
                              ">500", ">1k", ">5k", ">10k+"))+
  labs(x = "Sequence duplication levels", y = "% of library")

ggarrange(p1, p3, p5, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom", labels = c("B", "C", "D"))
ggarrange(p2, p4, p6, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom", labels = c("E", "F", "G"))


