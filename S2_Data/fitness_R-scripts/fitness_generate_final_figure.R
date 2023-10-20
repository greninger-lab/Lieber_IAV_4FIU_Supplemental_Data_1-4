library(tidyverse)
library(cowplot)
library(reshape2)
library(readxl)

setwd('/Volumes/lizso_backup_drive/Lieber_Flu_final/Figure-2F_Fitness/')
af_full <- read_csv("fitness_RAVA_output_to_figure_input/fitness_summary_lineage_mutations.csv")

samples <- as.data.frame(str_split_fixed(af_full$Passage, "-", 2))
colnames(samples) <- c("Virus", "PassageNum")

af_only <- af_full[,grepl("AF_", colnames(af_full))]

af <- cbind(samples, af_only)
af <- melt(af)
af$variable <- gsub("AF_", "", af$variable)
split <- as.data.frame(str_split_fixed(af$variable, ":", 2))
af <- cbind(af %>%select("Virus", "PassageNum"), split, af %>% select(value))
colnames(af)[3:5] <- c("Gene", "Mutation", "Allele Frequency")

af$Virus <- paste0("Virus ", af$Virus)
af[is.na(af)] <- 0

af$PassageNum <- gsub("Inoculum", "Inoc", af$PassageNum)

af$virus_mut <- paste(af$Virus, af$Gene, af$Mutation)
af$virus_mut <- gsub(" ", "_", af$virus_mut)
af$`Gene (Mutation)` <- paste0(af$Gene, " (", af$Mutation, ")")
afsum0 <- af
########################################################################################################################################
########################################################################################################################################
# MEDIAN IS GREATER THAN 0
########################################################################################################################################
########################################################################################################################################

muts <- af %>% group_by(virus_mut) %>% summarize(median=(median(`Allele Frequency`)))
write_csv(muts, "fitness_generate_final_figure/fitness_IAV_muts_all.csv")
muts <- muts %>% filter(median > 0) 

write_csv(af, "fitness_generate_final_figure/fitness_IAV_af.csv")
write_csv(af_full, "fitness_generate_final_figure/fitness_IAV_af_full.csv")
write_csv(muts, "fitness_generate_final_figure/fitness_IAV_muts_filteredMedian0.csv")

af <- af %>% filter(af$virus_mut %in% muts$virus_mut, PassageNum != "inoc")
af <- af %>% mutate(PassageNum = str_replace_all(PassageNum, "P",""),
                    `Allele Frequency` = round(`Allele Frequency`, digits = 2))
af$PassageNum <- factor(af$PassageNum, levels=c('0', '1', '2', '3', '4', '5'))


#superscript https://stackoverflow.com/questions/28978011/how-to-subscript-the-x-axis-tick-label
#fix legend spacing https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2

color_codes <- c("PA (N222S)"  = "cornflowerblue",
                 "PB2 (E180K)" = "red2",
                 "PB1 (M339I)" = "black",
                 "PB2 (Y488C)" = "red2",
                 "PB1 (V285I)" = "black",
                 "PB1 (M290V)" = "black",
                 "PA (M579I)"  = "cornflowerblue",
                 "PB1 (T46A)" = "black",
                 "PB2 (K189R)" = "red2" ,
                 "PB2 (E191K)" = "orange",
                 "PB2 (T491M)" = "orange",
                 "PA (S395N)" = "cornflowerblue", 
                 "   " = NA, 
                 "  " = NA, 
                 " " = NA)
shape_codes <- c("PA (N222S)"  = 1,
                 "PB2 (E180K)" = 2,
                 "PB1 (M339I)" = 3,
                 "PB2 (Y488C)" = 4,
                 "PB1 (V285I)" = 5,
                 "PB1 (M290V)" = 6,
                 "PA (M579I)"  = 7,
                 "PB1 (T46A)" = 8,
                 "PB2 (K189R)" = 9 ,
                 "PB2 (E191K)" = 10,
                 "PB2 (T491M)" = 11,
                 "PA (S395N)" = 12,
                 "   " = NA, 
                 "  " = NA, 
                 " " = NA)
#virus 1
df1 <- af %>% filter(Virus == "Virus lin1", PassageNum != "inoc")
legend_labels1 <- paste0(df1$`Gene (Mutation)`[df1$PassageNum==df1$PassageNum[max(as.numeric(df1$PassageNum))]],", ", df1$`Allele Frequency`[df1$PassageNum==df1$PassageNum[max(as.numeric(df1$PassageNum))]], "%")
df1$PassageNum <- factor(df1$PassageNum, levels=levels(af$PassageNum))
df1$`Gene (Mutation)` <- factor(df1$`Gene (Mutation)`, levels=c(unique(df1$`Gene (Mutation)`), "   ", "  ", " "))
A <- ggplot(df1, aes(x=PassageNum, y=`Allele Frequency`, color=`Gene (Mutation)`, group=`Gene (Mutation)`, shape = `Gene (Mutation)`)) + 
  geom_point(alpha=0.5, size=1.5) + geom_line() +
  scale_color_manual(values=color_codes, labels = c(legend_labels1, "", "", ""), na.value="white", drop=FALSE, name="") +
  scale_shape_manual(values=shape_codes, labels = c(legend_labels1, "", "", ""), na.value=1, drop=FALSE, name="")+
  labs(x="Passage number") +
  ylab("Allele Frequency (%)") +
  ggtitle("#1 - fitness test") +
  ylim(c(-1,105)) +
  theme_classic(base_size=6)+
  theme(legend.position="bottom", 
        legend.direction = "vertical",
        legend.spacing.x = unit(0, 'char'),
        axis.text.x = element_text(size=6),
        axis.title.x = element_text(size=4),
        legend.title = element_blank(),
        legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt")))
A


#virus 2
df2 <- af %>% filter(Virus == "Virus lin2", PassageNum != "inoc")
legend_labels2 <- paste0(df2$`Gene (Mutation)`[df2$PassageNum==df2$PassageNum[max(as.numeric(df2$PassageNum))]],", ", df2$`Allele Frequency`[df2$PassageNum==df2$PassageNum[max(as.numeric(df2$PassageNum))]], "%")
df2$PassageNum <- factor(df2$PassageNum, levels=levels(af$PassageNum))
df2$`Gene (Mutation)` <- factor(df2$`Gene (Mutation)`, levels=c(unique(df2$`Gene (Mutation)`), " "))
B <- ggplot(df2, aes(x=PassageNum, y=`Allele Frequency`, color=`Gene (Mutation)`, group=`Gene (Mutation)`, shape = `Gene (Mutation)`)) + 
  geom_point(alpha=0.5, size=1.5) + geom_line() +
  scale_color_manual(values=color_codes, labels = c(legend_labels2, ""), na.value="white", drop=FALSE, name="") +
  scale_shape_manual(values=shape_codes, labels = c(legend_labels2, "") , na.value=1, drop=FALSE, name="")+
  labs(x="Passage number") +
  ylab("Allele Frequency (%)") +
  ggtitle("#2 - fitness test") +
  ylim(c(-1,105)) +
  theme_classic(base_size=6)+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="bottom", 
        legend.direction = "vertical",
        legend.spacing.x = unit(0, 'char'),
        axis.text.x = element_text(size=6),
        axis.title.x = element_text(size=4),
        legend.title = element_blank(),
        legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt")))
B

#virus 3
df3 <- af %>% filter(Virus == "Virus lin3", PassageNum != "inoc")
legend_labels3 <- paste0(df3$`Gene (Mutation)`[df3$PassageNum==df3$PassageNum[max(as.numeric(df3$PassageNum))]],", ", df3$`Allele Frequency`[df3$PassageNum==df3$PassageNum[max(as.numeric(df3$PassageNum))]], "%")
df3$PassageNum <- factor(df3$PassageNum, levels=levels(af$PassageNum))
df3$`Gene (Mutation)` <- factor(df3$`Gene (Mutation)`, levels=c(unique(df3$`Gene (Mutation)`), " ", "  "))
C <- ggplot(df3, aes(x=PassageNum, y=`Allele Frequency`, color=`Gene (Mutation)`, group=`Gene (Mutation)`, shape = `Gene (Mutation)`)) + 
  geom_point(alpha=0.5, size=1.5) + geom_line() +
  scale_color_manual(values=color_codes, labels = c(legend_labels3, "", ""), na.value="white", drop=FALSE, name="") +
  scale_shape_manual(values=shape_codes, labels = c(legend_labels3, "", ""), na.value=1, drop=FALSE, name="")+
  labs(x="Passage number") +
  ylab("Allele Frequency (%)") +
  ggtitle("#3 - fitness test") +
  ylim(c(-1,105)) +
  theme_classic(base_size=6)+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="bottom", 
        legend.direction = "vertical",
        legend.spacing.x = unit(0, 'char'),
        axis.text.x = element_text(size=6),
        axis.title.x = element_text(size=4),
        legend.title = element_blank(),
        legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt")))
C

#virus 4
df4 <- af %>% filter(Virus == "Virus lin4", PassageNum != "inoc")
legend_labels4 <- paste0(df4$`Gene (Mutation)`[df4$PassageNum==df4$PassageNum[max(as.numeric(df4$PassageNum))]],", ", df4$`Allele Frequency`[df4$PassageNum==df4$PassageNum[max(as.numeric(df4$PassageNum))]], "%")
df4$PassageNum <- factor(df4$PassageNum, levels=levels(af$PassageNum))
df4$`Gene (Mutation)` <- factor(df4$`Gene (Mutation)`, levels=c(unique(df4$`Gene (Mutation)`), " "))
D <- ggplot(df4, aes(x=PassageNum, y=`Allele Frequency`, color=`Gene (Mutation)`, group=`Gene (Mutation)`, shape = `Gene (Mutation)`)) + 
  geom_point(alpha=0.5, size=1.5) + geom_line() +
  scale_color_manual(values=color_codes, labels = c(legend_labels4, ""), na.value="white", drop=FALSE, name="") +
  scale_shape_manual(values=shape_codes, labels = c(legend_labels4, ""), na.value=1, drop=FALSE, name="")+
  labs(x="Passage number") +
  ylab("Allele Frequency (%)") +
  ggtitle("#4 - fitness test") +
  ylim(c(-1,105)) +
  theme_classic(base_size=6)+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="bottom", 
        legend.direction = "vertical",
        legend.spacing.x = unit(0, 'char'),
        axis.text.x = element_text(size=6),
        axis.title.x = element_text(size=4),
        legend.title = element_blank(),
        legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt")))
D

#virus 5
df5 <- af %>% filter(Virus == "Virus lin5", PassageNum != "inoc")
legend_labels5 <- paste0(df5$`Gene (Mutation)`[df5$PassageNum==df5$PassageNum[max(as.numeric(df5$PassageNum))]],", ", df5$`Allele Frequency`[df5$PassageNum==df5$PassageNum[max(as.numeric(df5$PassageNum))]], "%")
df5$PassageNum <- factor(df5$PassageNum, levels=levels(af$PassageNum))
df5$`Gene (Mutation)` <- factor(df5$`Gene (Mutation)`, levels=c(unique(df5$`Gene (Mutation)`), " ", "  "))
E <- ggplot(df5, aes(x=PassageNum, y=`Allele Frequency`, color=`Gene (Mutation)`, group=`Gene (Mutation)`, shape = `Gene (Mutation)`)) + 
  geom_point(alpha=0.5, size=1.5) + geom_line() +
  scale_color_manual(values=color_codes, labels = c(legend_labels5, "", ""), na.value="white", drop=FALSE, name="") +
  scale_shape_manual(values=shape_codes, labels = c(legend_labels5, "", ""), na.value=1, drop=FALSE, name="")+
  labs(x="Passage number") +
  ylab("Allele Frequency (%)") +
  ggtitle("#5 - fitness test") +
  ylim(c(-1,105)) +
  theme_classic(base_size=6)+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="bottom", 
        legend.direction = "vertical",
        legend.spacing.x = unit(0, 'char'),
        axis.text.x = element_text(size=6),
        axis.title.x = element_text(size=4),
        legend.title = element_blank(),
        legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt")))
E

#virus 6
df6 <- af %>% filter(Virus == "Virus lin6", PassageNum != "inoc")
legend_labels6 <- paste0(df6$`Gene (Mutation)`[df6$PassageNum==df6$PassageNum[max(as.numeric(df6$PassageNum))]],", ", df6$`Allele Frequency`[df6$PassageNum==df6$PassageNum[max(as.numeric(df6$PassageNum))]], "%")
df6$PassageNum <- factor(df6$PassageNum, levels=levels(af$PassageNum))
df6$`Gene (Mutation)` <- factor(df6$`Gene (Mutation)`, levels=c(unique(df6$`Gene (Mutation)`), " "))

f <- ggplot(df6, aes(x=PassageNum, y=`Allele Frequency`, color=`Gene (Mutation)`, group=`Gene (Mutation)`, shape = `Gene (Mutation)`)) + 
  geom_point(alpha=0.5, size=1.5) + geom_line() +
  scale_color_manual(values=color_codes, labels = c(legend_labels6, ""), na.value="white", drop=FALSE, name="") +
  scale_shape_manual(values=shape_codes, labels = c(legend_labels6, ""), na.value=1, drop=FALSE, name="")+
  labs(x="Passage number") +
  ylab("Allele Frequency (%)") +
  ggtitle("#6 - fitness test") +
  ylim(c(-1,105)) +
  theme_classic(base_size=6) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="bottom", 
        legend.direction = "vertical",
        legend.spacing.x = unit(0, 'char'),
        axis.text.x = element_text(size=6),
        axis.title.x = element_text(size=4),
        legend.title = element_blank(),
        legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt")))
f

cow <- plot_grid(A, B, C, D, E, f, ncol=6, rel_widths = c(1.2,1,1,1,1,1))

pdf(file="fitness_generate_final_figure/Fig2F_fitness.pdf", width=6.5, height=2)
cow
dev.off()
