library(tidyverse)
library(cowplot)
library(reshape2)
library(readxl)

setwd('/Volumes/lizso_backup_drive/Lieber_Flu_final/Figure-4E_Mouse/')
af_full <- read_csv("mouse_RAVA_output_to_figure_input/mouse_summary_lineage_mutations.csv")
samples <- as.data.frame(str_split_fixed(af_full$Passage, "-", 3))
colnames(samples) <- c("Virus", "Input", "Animal")

af_only <- af_full[,grepl("AF_", colnames(af_full))]

af <- cbind(samples, af_only)
af <- melt(af)
af$variable <- gsub("AF_", "", af$variable)
split <- as.data.frame(str_split_fixed(af$variable, ":", 2))
af <- cbind(af %>%select("Virus", "Input", "Animal"), split, af %>% select(value))
colnames(af)[4:6] <- c("Gene", "Mutation", "Allele Frequency")

af$Input <- gsub("5E3", "5.0x10^3", af$Input)
af$Input <- gsub("1E3", "1.0x10^3", af$Input)
af$Input <- gsub("100", "1.0x10^2", af$Input)
af$Input <- gsub("1E2", "1.0x10^2", af$Input)
af$Input <- gsub("2.5E1", "2.5x10^1", af$Input)
af$Input <- gsub("25", "2.5x10^1", af$Input)


af$Virus <- paste0("Virus ", af$Virus)
af[is.na(af)] <- 0

af$Input <- gsub("inoc", "Inoc", af$Input)

af$virus_mut <- paste(af$Virus, af$Gene, af$Mutation)
af$virus_mut <- gsub(" ", "_", af$virus_mut)
af$`Gene (Mutation)` <- paste0(af$Gene, " (", af$Mutation, ")")

muts <- af %>% group_by(virus_mut) %>% summarize(median=(median(`Allele Frequency`)))
muts <- muts %>% filter(median > 0)

write_csv(muts, "mouse_generate_final_figure/mouse_IAV_muts_filter_med0.csv")
write_csv(af, "mouse_generate_final_figure/mouse_IAV_af.csv")
write_csv(af_full, "mouse_generate_final_figure/mouse_IAV_af_full.csv")

af <- af %>% filter(af$virus_mut %in% muts$virus_mut)

af$Input <- factor(af$Input, levels=c("Inoc", '5.0x10^3', '1.0x10^3', '5.0x10^2', '1.0x10^2', '2.5x10^1', '5'))


#superscript https://stackoverflow.com/questions/28978011/how-to-subscript-the-x-axis-tick-label
#fix legend spacing https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2

#virus 1
df1 <- af %>% filter(Virus == "Virus lin1")
df1$Input <- factor(df1$Input, levels=levels(af$Input))
df1$`Gene (Mutation)` <- factor(df1$`Gene (Mutation)`, levels=c(unique(df1$`Gene (Mutation)`), " ", ""))
A <- ggplot(df1, aes(x=Input, y=`Allele Frequency`, color=`Gene (Mutation)`)) + 
  geom_jitter(width=0.1, height=0.2, alpha=0.5) +
  scale_color_manual(values=c("orange", NA, NA), na.value="white", drop=FALSE, name="") +
  labs(x=expression(paste(TCID[50],  " units/mouse", sep=""))) +
  ylab("Allele Frequency (%)") +
  ggtitle("#1 - NGS lungs") +
  ylim(c(-1,101)) +
  theme_classic(base_size=6)+
  theme(legend.position="bottom", 
        legend.direction = "vertical",
        legend.spacing.x = unit(-0.5, 'char'),
        # legend.justification = c(0, 0),
        axis.text.x = element_text(size=6, angle = 45, hjust =1, vjust =1),
        legend.title = element_blank(),
        legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt")))
A


#virus 2
df2 <- af %>% filter(Virus == "Virus lin2")
df2$Input <- factor(df2$Input, levels=levels(af$Input))
df2$`Gene (Mutation)` <- factor(df2$`Gene (Mutation)`, levels=c(unique(df2$`Gene (Mutation)`)))
B <- ggplot(df2, aes(x=Input, y=`Allele Frequency`, color=`Gene (Mutation)`)) + 
  geom_jitter(width=0.1, height=0.2, alpha=0.5) +
  scale_color_manual(values=c("blue", "red3", "gold2"), name="") +
  labs(x=expression(paste(TCID[50],  " units/mouse", sep=""))) +
  ylab("Allele Frequency (%)") +
  ggtitle("#2 - NGS lungs") +
  ylim(c(-1,101)) +
  theme_classic(base_size=6)+
  theme(axis.title.y = element_blank(), 
       axis.text.y = element_blank(),
       axis.text.x = element_text(size=6, angle = 45, hjust =1,vjust=1),
        legend.position="bottom", 
       legend.direction="vertical",
        legend.spacing.x = unit(-0.5, 'char'),
       #legend.spacing.y = unit(-1, "char"),
       # legend.justification = c(0, 0),
        legend.title = element_blank(),
       legend.margin = margin(c(-10, 0, 0, 0)),
       legend.text = element_text(margin = margin(r = 5, unit = "pt"))) 
B

#virus 3
df3 <- af %>% filter(Virus == "Virus lin3")
df3$Input <- factor(df3$Input, levels=levels(af$Input))
df3$`Gene (Mutation)` <- factor(df3$`Gene (Mutation)`, levels=c(unique(df3$`Gene (Mutation)`), ""))
C <- ggplot(df3, aes(x=Input, y=`Allele Frequency`, color=`Gene (Mutation)`)) + 
  geom_jitter(width=0.1, height=0.2, alpha=0.5) +
  scale_color_manual(values=c("springgreen2", "orangered2", NA), na.value="white", drop=FALSE, name="") +
  labs(x=expression(paste(TCID[50],  " units/mouse", sep=""))) +
  ylab("Allele Frequency (%)") +
  ggtitle("#3 - NGS lungs") +
  ylim(c(-1,101)) +
  theme_classic(base_size=6)+
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=6, angle = 45, hjust =1,vjust=1),
        legend.position="bottom", 
        legend.direction="vertical",
        legend.spacing.x = unit(-0.5, 'char'),
        #legend.spacing.y = unit(-1, "char"),
        # legend.justification = c(0, 0),
          legend.title = element_blank(),
        legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt"))) 
C

#virus 4
df4 <- af %>% filter(Virus == "Virus lin4")
df4$Input <- factor(df4$Input, levels=levels(af$Input))
df4$`Gene (Mutation)` <- factor(df4$`Gene (Mutation)`, levels=c(unique(df4$`Gene (Mutation)`)))
D <- ggplot(df4, aes(x=Input, y=`Allele Frequency`, color=`Gene (Mutation)`)) + 
  geom_jitter(width=0.1, height=0.2, alpha=0.5) +
  scale_color_manual(values=c("magenta3", "darkturquoise", "darkorange2"), name="") +
  labs(x=expression(paste(TCID[50],  " units/mouse", sep=""))) +
  ylab("Allele Frequency (%)") +
  ggtitle("#4 - NGS lungs") +
  ylim(c(-1,101)) +
  theme_classic(base_size=6)+
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=6, angle = 45, hjust =1,vjust=1),
        legend.position="bottom", 
        legend.direction="vertical",
        legend.spacing.x = unit(-0.5, 'char'),
        #legend.spacing.y = unit(-1, "char"),
        # legend.justification = c(0, 0),
          legend.title = element_blank(),
        legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt"))) 
D

#virus 5
df5 <- af %>% filter(Virus == "Virus lin5")
df5$Input <- factor(df5$Input, levels=levels(af$Input))
df5$`Gene (Mutation)` <- factor(df5$`Gene (Mutation)`, levels=c(unique(df5$`Gene (Mutation)`), ""))
E <- ggplot(df5, aes(x=Input, y=`Allele Frequency`, color=`Gene (Mutation)`)) + 
  geom_jitter(width=0.1, height=0.2, alpha=0.5) +
  scale_color_manual(values=c("red", "orange", NA), na.value="white", drop=FALSE, name="") +
  labs(x=expression(paste(TCID[50],  " units/mouse", sep=""))) +
  ylab("Allele Frequency (%)") +
  ggtitle("#5 - NGS lungs") +
  ylim(c(-1,101)) +
  theme_classic(base_size=6)+
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=6, angle = 45, hjust =1,vjust=1),
        legend.position="bottom", 
        legend.direction="vertical",
        legend.spacing.x = unit(-0.5, 'char'),
        #legend.spacing.y = unit(-1, "char"),
        # legend.justification = c(0, 0),
          legend.title = element_blank(),
         legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt"))) 
E

#virus 6
df6 <- af %>% filter(Virus == "Virus lin6")
df6$Input <- factor(df6$Input, levels=levels(af$Input))
df6$`Gene (Mutation)` <- factor(df6$`Gene (Mutation)`, levels=c(unique(df6$`Gene (Mutation)`)))

f <- ggplot(df6, aes(x=Input, y=`Allele Frequency`, color=`Gene (Mutation)`)) + 
  geom_jitter(width=0.1, height=0.2, alpha=0.5) +
  scale_color_manual(values=c("chocolate3", "darkorchid3", "darkturquoise"), name="") +
  labs(x=expression(paste(TCID[50],  " units/mouse", sep=""))) +
  ylab("Allele Frequency (%)") +
  ggtitle("#6 - NGS lungs") +
  ylim(c(-1,101)) +
  theme_classic(base_size=6) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=6, angle = 45, hjust =1,vjust=1),
        legend.position="bottom", 
        legend.direction="vertical",
        legend.spacing.x = unit(-0.5, 'char'),
        #legend.spacing.y = unit(-1, "char"),
        # legend.justification = c(0, 0),
        legend.title = element_blank(),
        legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt"))) 
f

cow <- plot_grid(A, B, C, D, E, f, ncol=6, rel_widths = c(1.2,1,1,1,1,1))
#cow

#cow2 <- plot_grid(cow, legend, ncol=2, rel_widths=c(5, 1))


pdf(file="mouse_generate_final_figure/Figure4E_mouse.pdf", width=6.5, height=2.2)
cow
dev.off()








