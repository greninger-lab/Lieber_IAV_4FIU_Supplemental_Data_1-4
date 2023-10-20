library(tidyverse)
library(cowplot)
library(reshape2)
library(readxl)
library(dplyr)

setwd('/Volumes/lizso_backup_drive/Lieber_Flu_final/Figure-5H_Ferret-Transmission/')
af_full <- read_csv("ferret_RAVA_output_to_figure_input/ferret_summary_lineage_mutations.csv")

samples <- as.data.frame(str_split_fixed(af_full$Passage, "_", 3))
colnames(samples) <- c("Virus", "Input", "Animal")
samples$Virus <- str_replace_all(samples$Virus, "lin", "")
# samples$Input <- str_replace_all(samples$Input, )
af_only <- af_full[,grepl("AF_", colnames(af_full))]

af <- cbind(samples, af_only)
af <- melt(af)
af$variable <- gsub("AF_", "", af$variable)
split <- as.data.frame(str_split_fixed(af$variable, ":", 2))
af <- cbind(af %>%select("Virus", "Input", "Animal"), split, af %>% select(value))
colnames(af)[4:6] <- c("Gene", "Mutation", "Allele Frequency")

af$Virus <- paste0("Virus ", af$Virus)
af[is.na(af)] <- 0

af$virus_mut <- paste(af$Virus, af$Gene, af$Mutation)
af$virus_mut <- gsub(" ", "_", af$virus_mut)
af$`Gene (Mutation)` <- paste0(af$Gene, " (", af$Mutation, ")")

# filter mutations to plot for each lineage by sum of allele frequencies across ferrets
muts <- af %>% group_by(virus_mut) %>% summarize(sum=(sum(`Allele Frequency`)))


#superscript https://stackoverflow.com/questions/28978011/how-to-subscript-the-x-axis-tick-label
#fix legend spacing https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2

# include mutations where sum(allele frequency) across ferrets > 3%
muts <- muts %>% filter(sum > 3) 

write_csv(muts, "ferret_generate_final_figure/ferret_IAV_muts_filter_sum3.csv")

af <- af %>% filter(af$virus_mut %in% muts$virus_mut)

plot_af <- af
plot_af <- plot_af %>% add_row(Virus = "Virus 4", Input = 'contact-treated', `Gene (Mutation)` = "PA (S395N)")
plot_af <- plot_af %>% add_row(Virus = "Virus 5", Input = 'source-treated', `Gene (Mutation)` = "PA (N222S)")
plot_af <- plot_af %>% add_row(Virus = "Virus 5", Input = 'contact-vehicle', `Gene (Mutation)` = "PA (N222S)")
plot_af <- plot_af %>% add_row(Virus = "Virus 5", Input = 'contact-treated', `Gene (Mutation)` = "PA (N222S)")

plot_af$Input <- factor(plot_af$Input, levels=c("inoc", 'source-vehicle', 'source-treated', 'contact-vehicle', 'contact-treated'))

#virus 4
df4 <- plot_af %>% filter(Virus == "Virus 4")
df4$Input <- factor(df4$Input, levels=c(levels(plot_af$Input)))
df4$`Gene (Mutation)` <- factor(df4$`Gene (Mutation)`, levels=c(unique(df4$`Gene (Mutation)`), ""))
C <- ggplot(df4, aes(x=Input, y=`Allele Frequency`, color=`Gene (Mutation)`)) + 
  geom_jitter(width=0.1, height=0.2, alpha=0.5) +
  scale_color_manual(values=c("magenta3", "darkturquoise", "darkorange2", NA), na.value="white", drop=FALSE, name="") +
  labs(x="Treatment Group") +
  ylab("Allele Frequency (%)") +
  ggtitle("#4 - Ferret NGS") +
  # scale_x_discrete(labels=c("inoculum", "vehicle source", "treated source (12hpi)", "vehicle contact", "treated contact (12hpi)")) +
  ylim(c(-1,101)) +
  theme_classic(base_size=6)+
  theme(legend.position="bottom", 
        legend.direction = "horizontal",
        legend.spacing.x = unit(-0.5, 'char'),
        axis.text.x = element_text(size=5,angle=45, hjust=1),
        legend.title = element_blank(),
        legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt")))
C


#virus 5
df5 <- plot_af %>% filter(Virus == "Virus 5")
df5$Input <- factor(df5$Input, levels=levels(plot_af$Input))
df5$`Gene (Mutation)` <- factor(df5$`Gene (Mutation)`, levels=c(unique(df5$`Gene (Mutation)`), "", " "))
D <- ggplot(df5, aes(x=Input, y=`Allele Frequency`, color=`Gene (Mutation)`)) + 
  geom_jitter(width=0.1, height=0.2, alpha=0.5) +
  scale_color_manual(values=c("red", "orange", NA, NA), na.value="white", drop=FALSE, name="") +
  labs(x="Treatment Group") +
  ylab("Allele Frequency (%)") +
  ggtitle("#5 - Ferret NGS") +
  # scale_x_discrete(labels=c("inoculum", "vehicle source", "treated source (12hpi)", "vehicle contact", "treated contact (12hpi)")) +
  ylim(c(-1,101)) +
  theme_classic(base_size=6)+
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=5, angle=45, hjust=1),
        legend.position="bottom", 
        legend.direction="horizontal",
        legend.spacing.x = unit(-0.5, 'char'),
        legend.title = element_blank(),
        legend.margin = margin(c(-10, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 5, unit = "pt"))) 
D
cow <- plot_grid(C, D, ncol=2, rel_widths = c(1.05,1))

pdf(file="ferret_generate_final_figure/Figure5H_ferret-transmission.pdf", width=5.5, height=2.5)
cow
dev.off()

