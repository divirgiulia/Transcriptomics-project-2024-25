setwd("/Users/giuliadivirgilio/Desktop/BCG/I anno/II semestre/GENOMICS AND TRANSCRIPTOMICS/Transcriptomics/Pavesi_proj")
getwd
library(readxl)
library(tidyverse)
library(ggplot2)
library(recount3)
library(edgeR)
library(limma)
library(recount)
library(reshape2)


# loading the data (RDS format)
brain <- readRDS('rse_brain.RDS')
pancreas <- readRDS('rse_pancreas.RDS')
colon <- readRDS('rse_colon.RDS')

dim(brain) # 54042 genes, 2931 samples
dim(pancreas) # more than 54k genes, 360 samples
dim(colon) # 822 samples

# converting the mapping into counts, if the read is 50bp long and overlaps the gene, it counts as 50
assays(brain)$counts <- transform_counts(brain)
assays(pancreas)$counts <- transform_counts(pancreas)
assays(colon)$counts <- transform_counts(colon)

brain
pancreas
colon

# gives us information about the columns
colData(brain) #gtex means that the samples have info coming from the gtex portal
names(colData(brain))

# 
head(colData(brain)$gtex.smrin)

# n of annotations for each part of the sample
table(colData(brain)$gtex.smtsd) # how many times I find each of these annotations

# how many samples for each sex
table(colData(brain)$gtex.sex) # 2095 of sex1 and 836 of sex2

# age ranges
table(colData(brain)$gtex.age) # for privacy reasons

# number of input pair-end reads
head(colData(brain)$"recount_qc.star.number_of_input_reads_both")

head(colData(brain)$'recount_qc.star.uniquely_mapped_reads_%_both')

inputreads <- colData(brain)$"recount_qc.star.number_of_input_reads_both"
boxplot(inputreads)

minlib <- min(inputreads)
maxlib <- max(inputreads)
minlib
maxlib

# lots of outliers, there is some kind of problem - there was no quality control
mapped <- colData(brain)$"recount_qc.star.uniquely_mapped_reads_%_both"
boxplot(mapped)
min(mapped) # the lowest has 6.6% of pair end reads mapped

# 
min_map <- which.min(mapped)
colData(brain)[min_map,]

# RNA integrity number (should be >= 6)
rinplot <- colData(brain)$gtex.smrin
boxplot(rinplot)


head(colData(brain)$gtex.smrrnart)
boxplot(colData(brain)$gtex.smrrnart)

# % of mtRNA
head(colData(brain)$"recount_qc.aligned_reads%.chrm")  # it is better to exclude genes that are mapped from the mt
boxplot(colData(brain)$"recount_qc.aligned_reads%.chrm")
plot(colData(brain)$gtex.smrrnart, colData(brain)$"recount_qc.aligned_reads%.chrm")
# samples that have too much rRNA will also have too much mtRNA

# 25 diff piece of info associated with each gene
names(rowData(brain))
table(rowData(brain)$gbkey)
# Gene is a pseudoGene

rowRanges(brain)
# J_segment refers to the antibody

table(rowRanges(brain)@seqnames)

############### BRAIN
# selecting col 52, 54 and 55
colData(brain)$gtex.smrin[52] # 7.4
colData(brain)$gtex.smrin[53] # RIN < 6
colData(brain)$gtex.smrin[54] # 7.6
colData(brain)$gtex.smrin[55] # 7.5

# how many reads in each rep, in this case all above 30mln
colData(brain)$"recount_qc.star.number_of_input_reads_both"[52] #65043825
colData(brain)$"recount_qc.star.number_of_input_reads_both"[53] #43197924 x
colData(brain)$"recount_qc.star.number_of_input_reads_both"[54] #52130027
colData(brain)$"recount_qc.star.number_of_input_reads_both"[55] #35835118

# estimated fraction of rRNA, all < 0.1
colData(brain)$gtex.smrrnart[52] # 0.0345123
colData(brain)$gtex.smrrnart[53]
colData(brain)$gtex.smrrnart[54] # 0.00848132
colData(brain)$gtex.smrrnart[55] # 0.0786605

# mapped reads %, all above 80%
colData(brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[52] # 88.9
colData(brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[53] # 87
colData(brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[54] # 86.5
colData(brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[55] # 92

brain_selected <- brain[,c(52,54,55)]
counts_brain_selected <- assays(brain_selected)$counts


####### PANCREAS ########
# selecting col 52, 54 and 55
colData(pancreas)$gtex.smrin[52] # 6.4
colData(pancreas)$gtex.smrin[53] # 6.8
colData(pancreas)$gtex.smrin[54] # 6.5

# how many reads in each rep, in this case all above 30mln
colData(pancreas)$"recount_qc.star.number_of_input_reads_both"[52] # 33566561
colData(pancreas)$"recount_qc.star.number_of_input_reads_both"[53] # 43463677
colData(pancreas)$"recount_qc.star.number_of_input_reads_both"[54] # 40616305

# estimated fraction of rRNA, all < 0.1
colData(pancreas)$gtex.smrrnart[52] # 0.00470862
colData(pancreas)$gtex.smrrnart[53] # 0.00314337
colData(pancreas)$gtex.smrrnart[54] # 0.00326493

# mapped reads %, all above 80%
colData(pancreas)$"recount_qc.star.uniquely_mapped_reads_%_both"[52] # 83.8
colData(pancreas)$"recount_qc.star.uniquely_mapped_reads_%_both"[53] # 88.1
colData(pancreas)$"recount_qc.star.uniquely_mapped_reads_%_both"[54] # 90.4

pancreas_selected <- pancreas[,c(52,53,54)]
counts_pancreas_selected <- assays(pancreas_selected)$counts

####### COLON ########
# selecting col 52, 54 and 55
colData(colon)$gtex.smrin[52] # 6.9
colData(colon)$gtex.smrin[53]
colData(colon)$gtex.smrin[54] # 6.3
colData(colon)$gtex.smrin[55] # 7.4

# how many reads in each rep, in this case all above 30mln
colData(colon)$"recount_qc.star.number_of_input_reads_both"[52] # 34639257
colData(colon)$"recount_qc.star.number_of_input_reads_both"[54] # 37566922
colData(colon)$"recount_qc.star.number_of_input_reads_both"[55] # 38031370


# estimated fraction of rRNA, all < 0.1
colData(colon)$gtex.smrrnart[52] # 0.0191248
colData(colon)$gtex.smrrnart[53]
colData(colon)$gtex.smrrnart[54] # 0.00279837
colData(colon)$gtex.smrrnart[55] # 0.006781

# mapped reads %, all above 80%
colData(colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[52] # 88.5
colData(colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[53]
colData(colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[54] # 89.3
colData(colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[55] # 90.5


colon_selected <- colon[,c(52,54,55)]
counts_colon_selected <- assays(colon_selected)$counts

# the final count table with all three can be generated by joining them with cbind, e.g. final_count_table <- cbind(counts_tissue1_selected, counts_tissue2_selected, counts_tissue3_selected)
tissues_count_table <- cbind(counts_brain_selected, counts_pancreas_selected, counts_colon_selected)
View(tissues_count_table)

rownames(counts_brain_selected) <- rowData(brain)$gene_name
rownames(counts_pancreas_selected) <- rowData(pancreas)$gene_name
rownames(counts_colon_selected) <- rowData(colon)$gene_name

nrow(tissues_count_table)

# Moving onto finding DE genes
colnames(tissues_count_table) <- c("Brain52","Brain54","Brain55","Pancreas52", "Pancreas53","Pancreas54","Colon52","Colon54","Colon55")
y <- DGEList(counts=tissues_count_table)

group <- as.factor(c("Brain","Brain","Brain","Pancreas","Pancreas","Pancreas","Colon","Colon","Colon"))
y$samples$group <- group

y$samples$rin <- as.factor(c(colData(brain_selected)$gtex.smrin,colData(pancreas_selected)$gtex.smrin,colData(colon_selected)$gtex.smrin))
y$samples$slice <- as.factor(c(colData(brain_selected)$gtex.smtsd,colData(pancreas_selected)$gtex.smtsd,colData(colon_selected)$gtex.smtsd))
y$samples$sex <- as.factor(c(colData(brain_selected)$gtex.sex,colData(pancreas_selected)$gtex.sex,colData(colon_selected)$gtex.sex))
y$samples$age <- as.factor(c(colData(brain_selected)$gtex.age,colData(pancreas_selected)$gtex.age,colData(colon_selected)$gtex.age))
y$samples$rRNA <- as.factor(c(colData(brain_selected)$gtex.smrrnart,colData(pancreas_selected)$gtex.smrrnart,colData(colon_selected)$gtex.smrrnart))
y$samples$mapped <- as.factor(c(colData(brain_selected)$"recount_qc.star.uniquely_mapped_reads_%_both", colData(pancreas_selected)$"recount_qc.star.uniquely_mapped_reads_%_both",colData(colon_selected)$"recount_qc.star.uniquely_mapped_reads_%_both"))
y$samples$chrm <- as.factor(c(colData(brain_selected)$"recount_qc.aligned_reads%.chrm", colData(pancreas_selected)$"recount_qc.aligned_reads%.chrm",colData(colon_selected)$"recount_qc.aligned_reads%.chrm"))
y


table(rowSums(y$counts==0)==9)
# FALSE  TRUE 
# 39626 14416 

keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y) # we're left with more than 20k genes

logcpm_before <- cpm(y, log=TRUE) # we're storing in a vector the log of the counts per million before normalization
y <- calcNormFactors(y, method = "TMM")
y
logcpm_after <- cpm(y, log=TRUE)
boxplot(logcpm_before, notch=T, las=2)
boxplot(logcpm_after, notch=T, las=2)

#### cuter with ggplot
df_before <- melt(logcpm_before)
colnames(df_before) <- c("Gene", "Sample", "logCPM")
df_before$Group <- gsub("[0-9]+", "", df_before$Sample)
df_before$Group <- factor(df_before$Group, levels = c("Brain", "Pancreas", "Colon"))

df_after <- melt(logcpm_after)
colnames(df_after) <- c("Gene", "Sample", "logCPM")
df_after$Group <- gsub("[0-9]+", "", df_after$Sample)
df_after$Group <- factor(df_after$Group, levels = c("Brain", "Pancreas", "Colon"))

colors <- c("Brain" = "#FF0000", "Pancreas" = "#F2AD00", "Colon" = "#027ab0")
ggplot(df_before, aes(x = Sample, y = logCPM, fill = Group)) +
  geom_boxplot(notch = TRUE) +
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 12) +
  labs(title = "Log-CPM Before TMM Normalization",
       x = "Sample", y = "Log Counts Per Million") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# --- Boxplot for BEFORE normalization ---
ggplot(df_before, aes(x = Sample, y = logCPM, fill = Group)) +
  geom_boxplot(notch = TRUE) +
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 3)) +
  theme_minimal(base_size = 12) +
  labs(title = "Log-CPM Before TMM Normalization",
       x = "Sample", y = "Log Counts Per Million") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# --- Boxplot for AFTER normalization ---
ggplot(df_after, aes(x = Sample, y = logCPM)) +
  geom_boxplot(notch = TRUE, fill = "#238443") +
  theme_minimal(base_size = 12) +
  labs(title = "Log-CPM After TMM Normalization",
       x = "Sample", y = "Log Counts Per Million") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(df_after, aes(x = Sample, y = logCPM, fill = Group)) +
  geom_boxplot(notch = TRUE) +
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 12) +
  labs(title = "Log-CPM After TMM Normalization",
       x = "Sample", y = "Log Counts Per Million") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

##################

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design <- design[, c("Brain", "Pancreas", "Colon")]
design

logcpm <- cpm(y, log=TRUE)
plotMDS(logcpm, labels=group)
library(viridis)
library(wesanderson)
library(MetBrewer)

# Compute MDS coordinates manually
mds <- plotMDS(logcpm, plot = FALSE)
mds_df <- data.frame(
  Dim1 = mds$x,
  Dim2 = mds$y,
  Sample = colnames(logcpm),
  Group = group  # assuming 'group' is a vector or factor with same length as samples
)
ggplot(mds_df, aes(x = Dim1, y = Dim2, label = Sample, color = Group)) +
  geom_text(vjust = -0.5, size = 5, fontface = "bold") +
  theme_minimal(base_size = 12) +
  labs(
    title = "MDS Plot of Brain, Pancreas and Colon Samples (logCPM)",
    x = paste("Dimension 1 (", round(mds$var.explained[1], 1), "%)", sep = ""),
    y = paste("Dimension 2 (", round(mds$var.explained[2], 1), "%)", sep = "")
  ) +
  theme(legend.position = "right") +
  expand_limits(
    x = range(mds_df$Dim1) + c(-1, 1),
    y = range(mds_df$Dim2) + c(-1, 1)
  ) + scale_color_manual(values = wes_palette("Darjeeling1"))+  
  theme(plot.margin = margin(10, 30, 10, 30))
########################################################


mds_df2 <- data.frame(
  Dim1 = mds$x,
  Dim2 = mds$y,
  rRNA = as.character(y$samples$rRNA),      # labels
  chrm = as.character(y$samples$chrm),
  age = as.character(y$samples$age),
  sex = as.character(y$samples$sex),
  group = as.factor(y$samples$group)        # color by this
)

# Plot with color by group, label by rRNA

ggplot(mds_df2, aes(x = Dim1, y = Dim2, label = rRNA, color = group)) +
  geom_text(size = 4, fontface = "bold") +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  theme_minimal() +
  labs(color = "Sample Group", title = "MDS Plot of Sample Groups (by estimated % of rRNA reads)") +
  coord_cartesian(
    xlim = range(mds_df$Dim1) + c(-1, 1),  # expand x limits by 1 unit both sides
    ylim = range(mds_df$Dim2) + c(-1, 1)   # expand y limits by 1 unit both sides
  ) +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),   # top, right, bottom, left margins
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Plot with chrm labels
ggplot(mds_df2, aes(x = Dim1, y = Dim2, label = chrm, color = group)) +
  geom_text(size = 4, fontface = "bold") +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  theme_minimal() +
  labs(color = "Sample Group", title = "MDS Plot of Sample Groups (by estimated % of mtRNA reads)") +
  coord_cartesian(
    xlim = range(mds_df$Dim1) + c(-1, 1),  # expand x limits by 1 unit both sides
    ylim = range(mds_df$Dim2) + c(-1, 1)   # expand y limits by 1 unit both sides
  ) +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),   # top, right, bottom, left margins
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Plot with age labels
ggplot(mds_df2, aes(x = Dim1, y = Dim2, label = age, color = group)) +
  geom_text(size = 4, fontface = "bold") +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  theme_minimal() +
  labs(color = "Sample Group", title = "MDS Plot of Sample Groups (by age)") +
  coord_cartesian(
    xlim = range(mds_df$Dim1) + c(-1, 1),  # expand x limits by 1 unit both sides
    ylim = range(mds_df$Dim2) + c(-1, 1)   # expand y limits by 1 unit both sides
  ) +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),   # top, right, bottom, left margins
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Plot with sex labels
ggplot(mds_df2, aes(x = Dim1, y = Dim2, label = sex, color = group)) +
  geom_text(size = 4, fontface = "bold") +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  theme_minimal() +
  labs(color = "Sample Group", title = "MDS Plot of Sample Groups (by sex)") +
  coord_cartesian(
    xlim = range(mds_df$Dim1) + c(-1, 1),  # expand x limits by 1 unit both sides
    ylim = range(mds_df$Dim2) + c(-1, 1)   # expand y limits by 1 unit both sides
  ) +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),   # top, right, bottom, left margins
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

plotMDS(logcpm, labels=y$samples$rRNA)
plotMDS(logcpm, labels=y$samples$chrm) ##
plotMDS(logcpm, labels=y$samples$age)
plotMDS(logcpm, labels=y$samples$sex)  ##

#######

y <- estimateDisp(y, design)
plotBCV(y)
title("Biological Coefficient of Variation Plot")

# ggplot2 - BCV plot
#plot_bcv(y)
bcv_data <- data.frame(
  AveLogCPM = y$AveLogCPM,
  BCV = sqrt(y$tagwise.dispersion),
  Trend = sqrt(y$trended.dispersion)
)

common_bcv <- sqrt(y$common.dispersion)

ggplot(bcv_data, aes(x = AveLogCPM, y = BCV)) +
  geom_point(alpha = 0.4) +
  labs(
    title = "Biological Coefficient of Variation (BCV) Plot",
    x = "Average Log CPM",
    y = "Biological Coefficient of Variation"
  ) +  
  geom_point(alpha = 0.4, color = "black") + # BCV points
  geom_line(aes(y = Trend), color = "blue", size = 1) + # Trend line
  geom_hline(yintercept = common_bcv, color = "red") +
  theme_minimal()

fit <- glmQLFit(y, design)

# pancreas (top) vs brain (bottom)
qlfPB <- glmQLFTest(fit, contrast=c(-1,1,0))
#### brain vs pancreas ?
qlfBP <- glmQLFTest(fit, contrast=c(1,-1,0))
# colon (top) vs brain (bottom)
qlfCB <- glmQLFTest(fit, contrast=c(-1,0,1))
# colon (top) vs pancreas (bottom)
qlfCP <- glmQLFTest(fit, contrast=c(0,-1,1))


qlfPB
qlfCB
qlfCP

###### Pancreas vs everything else

qlfPC <- glmQLFTest(fit, contrast=c(0,1,-1))
qlfPC
summary(decideTests(qlfPC, p.value=0.05, adjust.method = "BH", lfc=1))
dtPC <- decideTests(qlfPC, p.value = 0.05, adjust.method = "BH", lfc = 1)
dtPB <- decideTests(qlfPB, p.value = 0.05, adjust.method = "BH", lfc = 1)
summary(dtPC)
summary(dtPB)
upPC <- dtPC[,1] == 1  # Pancreas vs Colon: up in Pancreas
upPB <- dtPB[,1] == 1  # Pancreas vs Brain: up in Pancreas
up_in_pancreas <- upPC & upPB
sum(up_in_pancreas)
up_gene_names <- rownames(dtPC)[up_in_pancreas]
head(up_gene_names)
up_gene_names
which(rowData(pancreas)$gene_name == "MFAP2")
boxplot(assays(brain)$TPM[6660,],
        assays(pancreas)$TPM[6660,], 
        assays(colon)$TPM[6660,], 
        outline=F, 
        names=c("Brain","Pancreas","Colon"))

# Brain vs everything else
qlfBP <- glmQLFTest(fit, contrast=c(1,-1,0))
qlfBC <- glmQLFTest(fit, contrast=c(1,0,-1))

summary(decideTests(qlfBP, p.value=0.05, adjust.method = "BH", lfc=1))
dtBP <- decideTests(qlfBP, p.value = 0.05, adjust.method = "BH", lfc = 1)
dtBC <- decideTests(qlfBC, p.value = 0.05, adjust.method = "BH", lfc = 1)
summary(dtBP)
summary(dtBC)
upBP <- dtBP[,1] == 1  # Brain vs Pancreas: up in Brain
upBC <- dtBC[,1] == 1  # Brain vs Colon: up in Brain
up_in_brain <- upBP & upBC
sum(up_in_brain) # 2471 genes up-reg in Brain
up_gene_names_brain <- rownames(dtBP)[up_in_brain]
head(up_gene_names_brain)
up_gene_names_brain
which(rowData(brain)$gene_name == "TMEM88B") # row 6235
boxplot(assays(brain)$TPM[6235,],
        assays(pancreas)$TPM[6235,], 
        assays(colon)$TPM[6235,], 
        outline=F, 
        names=c("Brain","Pancreas","Colon"))
wilcox.test(assays(brain)$TPM[6253,], assays(pancreas)$TPM[6253,])



# Colon vs everything else
qlfCB <- glmQLFTest(fit, contrast=c(-1,0,1))
qlfCP <- glmQLFTest(fit, contrast=c(0,-1,1))

dtCB <- decideTests(qlfCB, p.value = 0.05, adjust.method = "BH", lfc = 1)
dtCP <- decideTests(qlfCP, p.value = 0.05, adjust.method = "BH", lfc = 1)
summary(dtCB)
summary(dtCP)
upCB <- dtCB[,1] == 1  # Brain vs Pancreas: up in Brain
upCP <- dtCP[,1] == 1  # Brain vs Colon: up in Brain
up_in_colon <- upCB & upCP
sum(up_in_colon) # 431 genes up-reg in Brain
up_gene_names_colon <- rownames(dtCB)[up_in_colon]
head(up_gene_names_colon)
which(rowData(brain)$gene_name == "MXRA8") # row 6245
#

boxplot(assays(brain)$TPM[6245,],
        assays(pancreas)$TPM[6245,], 
        assays(colon)$TPM[6245,], 
        outline=F, 
        names=c("Brain","Pancreas","Colon"))
wilcox.test(assays(colon)$TPM[6245,], assays(brain)$TPM[6245,])
wilcox.test(assays(colon)$TPM[6245,], assays(pancreas)$TPM[6245,])

# Counts
num_up_pancreas <- sum(up_in_pancreas)
num_up_brain    <- sum(up_in_brain)
num_up_colon    <- sum(up_in_colon)

df <- data.frame(
  tissue   = c("Pancreas", "Brain", "Colon"),
  up_genes = c(num_up_pancreas, num_up_brain, num_up_colon)
)

library(ggplot2)

ggplot(df, aes(x = tissue, y = up_genes, fill = tissue)) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = up_genes),
    vjust = -0.4,
    size = 5
  ) +
  scale_fill_manual(values = colors, guide = "none") +
  labs(
    x = "Tissue",
    y = "Number of genes",
    title = "Number of Upregulated Genes per Tissue"
  ) +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, max(df$up_genes) * 1.1))



####### PANCREAS VS BRAIN
head(qlfPB$table)
topTags(qlfPB, n=10, adjust.method = "BH", sort.by = "PValue")
resultsPB <- topTags(qlfPB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)

# saving the results in a table
write.table(resultsPB, "resultsPB.txt")
summary(decideTests(qlfPB, p.value=0.05, adjust.method = "BH", lfc=1))

###### COLON VS BRAIN 
# 
head(qlfCB$table)
topTags(qlfCB, n=10, adjust.method = "BH", sort.by = "PValue")
resultsCB <- topTags(qlfCB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
summary(decideTests(qlfCB, p.value=0.05, adjust.method = "BH", lfc=1))
write.table(resultsCB, "resultsCB.txt")

#### COLON VS PANCREAS
head(qlfCP$table)
topTags(qlfCP, n=10, adjust.method = "BH", sort.by = "PValue")
resultsCP <- topTags(qlfCP, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
summary(decideTests(qlfCP, p.value=0.05, adjust.method = "BH", lfc=1))
write.table(resultsCP, "resultsCP.txt")


assays(brain)$TPM <- recount::getTPM(brain)
assays(colon)$TPM <- recount::getTPM(colon)
assays(pancreas)$TPM <- recount::getTPM(pancreas)
which(rowData(brain)$gene_name == "TMEM52")

boxplot(assays(brain)$TPM[6278,],
        assays(pancreas)$TPM[6278,], 
        assays(colon)$TPM[6278,], 
        outline=F, 
        names=c("Brain","Pancreas","Colon"))


library(dplyr)
library(tibble)
brain_tpm <- assays(brain)$TPM[6245, ]
pancreas_tpm <- assays(pancreas)$TPM[6245, ]
colon_tpm <- assays(colon)$TPM[6245, ]

df_MXRA8 <- tibble(
  TPM = c(brain_tpm, pancreas_tpm, colon_tpm),
  Tissue = factor(rep(c("Brain", "Pancreas", "Colon"), 
                      times = c(length(brain_tpm), length(pancreas_tpm), length(colon_tpm))),
                  levels = c("Brain", "Pancreas", "Colon"))  # <-- This sets the desired order
)

colors <- c("Brain" = "#FF0000", "Pancreas" = "#F2AD00", "Colon" = "#00A08A")
ggplot(df_MXRA8, aes(x = Tissue, y = TPM, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(0, 200)) +  # Adjust 200 to your desired upper limit
  theme_minimal() +
  labs(title = "Expression of Gene MXRA8 (TPM)",
       y = "TPM", x = NULL) +
  theme(legend.position = "none")

# comparing expression levels in brain vs pancreas and brain vs colon w/ non-parametric test
wilcox.test(assays(brain)$TPM[6235,], assays(pancreas)$TPM[6235,]) # p-value < 2.2e-16
wilcox.test(assays(colon)$TPM[6278,], assays(pancreas)$TPM[6278,]) # p-value < 2.2e-16
wilcox.test(assays(brain)$TPM[6235,], assays(colon)$TPM[6235]) # p-value < 2.2e-16

