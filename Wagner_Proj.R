##------------------------------------------Wagner Project----------------------------------------------------------

library(dplyr)
library(tidyr)
library(limma)
library(openxlsx)
library(readr)
library(edgeR)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(ggrepel)
library(viridis)
library(pheatmap)
library(grid)
library(VennDiagram)
library(fgsea)
library(gggsea)
library(purrr)
library(ggvenn)

##---------------------------------------Data cleaning------------------------------------------------------

## Load the sample data 
sample_info <- read.xlsx("FromWagner_sample_mapping_updated_usethis_111524.xlsx")
raw_counts <- read.delim("FromJarek_annotated_counts_080724.txt", sep = '\t')

study1_raw_counts <- raw_counts[, c(1,4, which(colnames(raw_counts) == "WD01"): which(colnames(raw_counts) == "WD22"))]
study2_raw_counts <- raw_counts[, c(1,4, which(colnames(raw_counts) == "WD23"):which(colnames(raw_counts) == "WD52"))]

any(study1_raw_counts[, 2] == "") # check if there any blank values in study1 data
any(study2_raw_counts[,2] == "")  # check if there are any blank values in study2 data

## replace the blank values found in gene_names column for both studies using the ensemble id values 
study1_raw_counts[, 2][study1_raw_counts[, 2] == ""] <- study1_raw_counts[, 1][study1_raw_counts[, 2] == ""]
study2_raw_counts[,2][study2_raw_counts[,2] == ""] <- study2_raw_counts[,1][study2_raw_counts[,2] == ""]

study1_sample_data <- sample_info[1:22,]
study2_sample_data <- sample_info[25:54,]
rownames(study2_sample_data) <- NULL

##---------------------------------------Main Analysis - QC ---------------------------------------------------

#check if there are any duplicated genes in raw count data 
any(duplicated(study1_raw_counts$gene_name)) # found 1624 genes

dup_genes <- study1_raw_counts[duplicated(study1_raw_counts$gene_name), ] # add the duplicated genes data into a new df
dup_genes <-  dup_genes %>%                          # concatenate the gene name with the ensembl_id
  mutate(gene_name = paste0(ensembl_gene_id, "_", gene_name))

#replace the modified duplicated gene name values in the raw counts data
study1_raw_counts$gene_name <- ifelse(
  rownames(study1_raw_counts) %in% rownames(dup_genes),  # Check if row names match
  dup_genes$gene_name[match(rownames(study1_raw_counts), rownames(dup_genes))],  # Replace with modified names
  study1_raw_counts$gene_name  # Keep original name if no match
)

## create a DGE list for study1 data
study1_raw_counts <- study1_raw_counts[,-1]
study1_DGE <- DGEList(counts = study1_raw_counts)

## TMM normalization
study1_DGE <- calcNormFactors(study1_DGE)

## log2 transformation of count data and prior count = 3
study1_logCPM <- cpm(study1_DGE$counts, log = T, prior.count = 3) # logcpm values without gene names

study1_logCPM_file <- data.frame(Gene = study1_DGE$genes$gene_name, study1_logCPM) #logcpm values with gene names

## add group data to DGE list
study1_groups <- study1_sample_data %>% select(c(1,5))
study1_groups$Raw.count.file <- sub("WD(\\d)$", "WD0\\1", study1_groups$Raw.count.file) # adjusting the raw.count.file values so that they match with the DGE list colnames
study1_groups$Groups <- ifelse(               # extract group names from sample id column
  grepl("TNF", study1_groups$SAMPLE.ID), "TNF",
  ifelse(grepl("TCZ", study1_groups$SAMPLE.ID), "TCZ", NA)
)
study1_groups$Time <- ifelse(               # extract timepoints from sample id column
  grepl("pre", study1_groups$SAMPLE.ID), "pre",
  ifelse(grepl("post", study1_groups$SAMPLE.ID), "post", NA)
)

study1_groups$Patient <- sub("(\\w+).*", "\\1", study1_groups$SAMPLE.ID)
S1_unique_samples <- unique(study1_groups$Patient) 
study1_groups$Patient <- paste0("sub", match(study1_groups$Patient, S1_unique_samples)) ## extract the patient data from the sample id

study1_groups$Patient_Time <- paste0("sub", match(study1_groups$Patient, unique(study1_groups$Patient)), "_", study1_groups$Time,  "_", study1_groups$Raw.count.file) #combine patient and time data together

study1_groups <- study1_groups %>% group_by(Patient) %>% mutate(Patient_Combined = paste0(Patient, "_", paste(unique(Raw.count.file), collapse = "_"))
  ) %>%
  ungroup() # combine the patient and raw count file data

# replace the group values in the DGE list with the actual group values from the sample data
study1_DGE$samples$group <- study1_groups$Groups

# Filter the lowly expressed genes
keep <- rowSums(cpm(study1_DGE) >1) >= 11
study1_DGE <- study1_DGE[keep,,keep.lib.sizes = F]
length(rownames(study1_DGE$counts)) #12395 genes remained after filtering


## Density plot
S1_L <- mean(study1_DGE$samples$lib.size) * 1e-6
S1_M <- median(study1_DGE$samples$lib.size) * 1e-6
S1_logCPM_thres <- log2(10/S1_M + 2/S1_L)

#jpeg("Study1_DensityPlot.jpeg", width = 2500, height = 1800, res = 300)
#Density Plot 
S1_samples <- ncol(study1_DGE)
col <- brewer.pal(S1_samples, "Dark2")
par(mfrow=c(1,2))
plot(density(study1_logCPM[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main=" Raw_data", xlab="Log-cpm")
abline(v= S1_logCPM_thres, lty=3)
for (i in 2:S1_samples){
  den <- density(study1_logCPM[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(study1_DGE), text.col=col, bty="n")
study1_logCPM <- cpm(study1_DGE, log=TRUE, prior.count = 3)
plot(density(study1_logCPM[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main=" Filtered_data", xlab="Log-cpm")
abline(v=S1_logCPM_thres, lty=3)
for (i in 2:S1_samples){
  den <- density(study1_logCPM[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(study1_DGE), text.col=col, bty="n")
#dev.off()


## create MDS plot to check the variability in gene expr data
#jpeg("Study1_MDSplot.jpeg", width = 2500, height = 1800, res = 300)
plotMDS(study1_logCPM, labels = study1_DGE$samples$group, col = ifelse(study1_DGE$samples$group == "TNF", "blue", "red"), cex = 1)
#dev.off()

## PCA plot 
## Principal Component Analysis 
S1_data <- data.table(study1_logCPM) # convert logcpm matrix to a data table obj
S1_tdata <- data.table::transpose(S1_data) # transpose the data table
colnames(S1_tdata) <- study1_DGE$genes$gene_name #set colnames of the transposed data
rownames(S1_tdata) <- colnames(study1_logCPM) #set rownames of the transposed data

S1_pca <- prcomp(S1_tdata, scale. = T) # pca on transposed data

##screeplot 
screeplot(S1_pca)

#dev.off()

## proportion of variance explained by each component
summary(S1_pca)

S1_pca_data <- data.frame(S1_pca$x) # convert principal component values to a dataframe obj
S1_pca_data$Group <- study1_DGE$samples$group #Add the groups column to pca dataframe
S1_pca_data$Patient_Time <- study1_groups$Patient_Time

# Assign colors to groups
S1_pca_groups <- unique(study1_DGE$samples$group)
S1_colors <- brewer.pal(length(S1_pca_groups), "Set1")
S1_pca_data$Color <- S1_colors[as.numeric(factor(S1_pca_data$Group, levels = S1_pca_groups))]

#jpeg("Study1_PCA.jpeg", width = 2500, height = 1800, res = 300)
# Plot PCA
ggplot(data= S1_pca_data, aes(x=PC1, y=PC2)) +
  geom_point(size=5, aes(color=Group, shape = S1_Time ))+
  geom_hline(yintercept=0, colour="black", linetype="dashed") +
  geom_vline(xintercept=0, colour="black", linetype="dashed") +
  geom_text_repel(label= S1_pca_data$Patient_Time, colour="black", size=3, nudge_x=0.1, nudge_y=0.1, max.overlaps = Inf) +  
  xlab("PC1 (28%)")+
  ylab("PC2(9%)")+ 
  theme(axis.text=element_text(size=14, colour="black"),axis.title=element_text(size=14, face="bold", colour="black")) + 
  theme(panel.background = element_rect(fill = "white", colour = 'black'))+
  ggtitle("Study1_PCA_plot_Wagner_Study")
#dev.off()

##Differential gene expression analysis - Paired data
S1_subject <- factor(study1_groups$Patient_Combined)
S1_subject <- factor(S1_subject, levels = unique(S1_subject)) #get the subject data

## Time 
S1_Time <- factor(study1_groups$Time)
S1_Time <- factor(study1_groups$Time, levels = unique(S1_Time))

S1_paired_data <- data.frame(Samplenames = colnames(study1_DGE$counts), Time = S1_Time, Subject = S1_subject)

## create a design matrix for the time comparison
S1_design <- model.matrix(~0+Time+Subject, data = S1_paired_data)
colnames(S1_design) <- gsub("Time", "", colnames(S1_design))
colnames(S1_design) <- gsub("Subject", "", colnames(S1_design))
rownames(S1_design) <- S1_paired_data$Samplenames

## contrast matrix for time comparison
S1_contr.mtrx <- makeContrasts(
  Time = post - pre,
  levels = S1_design)
S1_contr.mtrx

## run voom
S1_v <- voom(study1_DGE, S1_design, plot = T)
S1_fit <- lmFit(S1_v, S1_design)
S1_fit <- contrasts.fit(S1_fit, S1_contr.mtrx)
S1_efit <- eBayes(S1_fit)

PostvsPre_ALL <- topTable(S1_efit, n = Inf, adjust.method = "fdr", sort.by = "P") # no fdr threshold was given to this data - contains all genes
PostvsPre_FDR_0.05 <- topTable(S1_efit, n = Inf, adjust.method = "fdr", sort.by = "P", p.value = 0.05) #found 4814 genes at FDR 0.05
PostvsPre_FDR_0.01 <- topTable(S1_efit, n = Inf, adjust.method = "fdr", sort.by = "P", p.value = 0.01) #found 2967 genes at FDR 0.01


## Differential gene expression analysis - Unpaired data
Group.Time <- factor(paste(study1_groups$Groups, study1_groups$Time, sep = "_"))
Group.Time <- factor(Group.Time, levels = unique(Group.Time))
S1_group <- factor(study1_groups$Groups)
S1_group <- factor(S1_group, levels = unique(S1_group))

S1_unpaired_data <- data.frame(Samplenames = colnames(study1_DGE$counts), Time = S1_Time, Group = S1_group, Group.Time = Group.Time)


#design matrix for unpaired analysis
S1_design1 <- model.matrix(~0+Group.Time, data = S1_unpaired_data)
colnames(S1_design1) <- gsub("Group.Time", "", colnames(S1_design1))
rownames(S1_design1) <- S1_unpaired_data$Samplenames

# contrast matrix for unpaired analysis
S1_contr.mtrx1 <- makeContrasts(
  TNFvsTCZ = (TNF_pre + TNF_post)/2 - (TCZ_pre + TCZ_post)/2,
  TNFvsTCZ_pre = TNF_pre - TCZ_pre,
  TNFvsTCZ_post = TNF_post - TCZ_post,
  TNFvsTCZ_Intxn = (TNF_pre - TNF_post) - (TCZ_pre - TCZ_post),
  levels = S1_design1)

## Run voom
S1_v1 <- voom(study1_DGE, S1_design1, plot = T)
S1_fit1 <- lmFit(S1_v1, S1_design1)
S1_fit1 <- contrasts.fit(S1_fit1, S1_contr.mtrx1)
S1_efit1 <- eBayes(S1_fit1)

TNFvsTCZ_ALL <- topTable(S1_efit1, n = Inf, coef = 1, adjust.method = "fdr", sort.by = "P") # no fdr threshold was given to this data
TNFvsTCZ_pre <- topTable(S1_efit1, n = Inf, coef = 2, adjust.method = "fdr", sort.by = "P") # no fdr threshold was given to this data
TNFvsTCZ_post <- topTable(S1_efit1, n = Inf, coef = 3, adjust.method = "fdr", sort.by = "P") # no fdr threshold was given to this data
TNFvsTCZ_Intxn <- topTable(S1_efit1, n = Inf, coef = 4, adjust.method = "fdr", sort.by = "P") # no fdr threshold was given to this data

topTable(S1_efit1, n= Inf, coef = 1, adjust.method = "fdr", sort.by = "P", p.value = 0.05) # no genes were found at FDR 0.05 for TNFvsTCZ comparison
topTable(S1_efit1, n = Inf, coef = 2, adjust.method = "fdr", sort.by = "P", p.value = 0.05) # no genes were found at FDR 0.05 for TNFvsTCZ_pre comparison
topTable(S1_efit1, n = Inf, coef = 3, adjust.method = "fdr", sort.by = "P", p.value = 0.05) # no genes were found at FDR 0.05 for TNFvsTCZ_post comparison
topTable(S1_efit1, n = Inf, coef = 4, adjust.method = "fdr", sort.by = "P", p.value = 0.05) # no genes were found at FDR 0.05 for TNFvsTCZ_Intxn comparison

topTable(S1_efit1, n= Inf, coef = 1, adjust.method = "fdr", sort.by = "P", p.value = 0.1) # no genes were found at FDR 0.1 for TNFvsTCZ comparison
topTable(S1_efit1, n = Inf, coef = 2, adjust.method = "fdr", sort.by = "P", p.value = 0.1) # no genes were found at FDR 0.1 for TNFvsTCZ_pre comparison
topTable(S1_efit1, n = Inf, coef = 3, adjust.method = "fdr", sort.by = "P", p.value = 0.1) # no genes were found at FDR 0.1 for TNFvsTCZ_post comparison
topTable(S1_efit1, n = Inf, coef = 4, adjust.method = "fdr", sort.by = "P", p.value = 0.1) # no genes were found at FDR 0.1 for TNFvsTCZ_Intxn comparison

topTable(S1_efit1, n= Inf, coef = 1, adjust.method = "fdr", sort.by = "P", p.value = 0.2) # no genes were found at FDR 0.2 for TNFvsTCZ comparison
topTable(S1_efit1, n = Inf, coef = 2, adjust.method = "fdr", sort.by = "P", p.value = 0.2) # no genes were found at FDR 0.2 for TNFvsTCZ_pre comparison
topTable(S1_efit1, n = Inf, coef = 3, adjust.method = "fdr", sort.by = "P", p.value = 0.2) # no genes were found at FDR 0.2 for TNFvsTCZ_post comparison
topTable(S1_efit1, n = Inf, coef = 4, adjust.method = "fdr", sort.by = "P", p.value = 0.2) # no genes were found at FDR 0.2 for TNFvsTCZ_Intxn comparison



##-------------------------------------------Removal of Outliers-------------------------------------------------------

study1_raw_counts_out <- subset(study1_raw_counts, select = -c(WD03, WD04, WD21, WD22))
study1_DGE_out <- DGEList(counts = study1_raw_counts_out) # convert the new raw counts data into DGE list

##TMM normalization
study1_DGE_out <- calcNormFactors(study1_DGE_out)

## log2 transformation of count data
study1_logCPM_out <- cpm(study1_DGE_out$counts, log = T, prior.count = 3)

study1_logCPM_out_file <- data.frame(Gene = study1_DGE_out$genes$gene_name, study1_logCPM_out)
#write.xlsx(study1_logCPM_out_file, "study1_logCPM.xlsx")

study1_groups_out <- data.frame(Groups = study1_groups$Groups[match(colnames(study1_DGE_out), study1_groups$Raw.count.file)]) #pull the group values that match the existing sample data

## add the group values to DGE list
study1_DGE_out$samples$group <- study1_groups_out$Groups

## Filter the low quality genes
# Filter the lowly expressed genes
keep1 <- rowSums(cpm(study1_DGE_out) >1) >= 9
study1_DGE_out <- study1_DGE_out[keep1,,keep.lib.sizes = F]
length(rownames(study1_DGE_out$counts)) #12363 genes remained after filtering


## Density plot
S1_L1 <- mean(study1_DGE_out$samples$lib.size) * 1e-6
S1_M1 <- median(study1_DGE_out$samples$lib.size) * 1e-6
S1_logCPM_thres1 <- log2(10/S1_M1 + 2/S1_L1)

#jpeg("Study1_DensityPlot_without_outliers.jpeg", width = 2500, height = 1800, res = 300)
#Density Plot 
S1_samples1 <- ncol(study1_DGE_out)
col <- brewer.pal(S1_samples1, "Dark2")
par(mfrow=c(1,2))
plot(density(study1_logCPM_out[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main=" Raw_data", xlab="Log-cpm")
abline(v= S1_logCPM_thres1, lty=3)
for (i in 2:S1_samples1){
  den <- density(study1_logCPM_out[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(study1_DGE_out), text.col=col, bty="n")
study1_logCPM_out <- cpm(study1_DGE_out, log=TRUE, prior.count = 3)
plot(density(study1_logCPM_out[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main=" Filtered_data", xlab="Log-cpm")
abline(v=S1_logCPM_thres1, lty=3)
for (i in 2:S1_samples1){
  den <- density(study1_logCPM_out[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(study1_DGE_out), text.col=col, bty="n")
#dev.off()

## PCA plot 
## Principal Component Analysis 
S1_data1 <- data.table(study1_logCPM_out) # convert logcpm matrix to a data table obj
S1_tdata1 <- data.table::transpose(S1_data1) # transpose the data table
colnames(S1_tdata1) <- study1_DGE_out$genes$gene_name #set colnames of the transposed data
rownames(S1_tdata1) <- colnames(study1_logCPM_out) #set rownames of the transposed data

S1_pca_out <- prcomp(S1_tdata1, scale. = T) # pca on transposed data

##screeplot 
screeplot(S1_pca_out)

#dev.off()

## proportion of variance explained by each component
summary(S1_pca_out)

S1_pca_data1 <- data.frame(S1_pca_out$x) # convert principal component values to a dataframe obj
S1_pca_data1$Group <- study1_DGE_out$samples$group #Add the groups column to pca dataframe

Patient_Time_out <- study1_groups$Patient_Time[match(colnames(study1_DGE_out), study1_groups$Raw.count.file)] ## adjust the patient data according to outlier removal

S1_pca_data1$Patient_Time <- Patient_Time_out # add the updated patient data to pca data1

# Assign colors to groups
S1_pca_groups1 <- unique(study1_DGE_out$samples$group)
S1_colors1 <- brewer.pal(length(S1_pca_groups1), "Set1")
S1_pca_data1$Color <- S1_colors1[as.numeric(factor(S1_pca_data1$Group, levels = S1_pca_groups1))]

#jpeg("Study1_PCA_without_outliers.jpeg", width = 2500, height = 1800, res = 300)
# Plot PCA
ggplot(data= S1_pca_data1, aes(x=PC1, y=PC2)) +
  geom_point(size=5, aes(color=Group, shape = S1_Time_out ))+
  geom_hline(yintercept=0, colour="black", linetype="dashed") +
  geom_vline(xintercept=0, colour="black", linetype="dashed") +
  geom_text_repel(label= S1_pca_data1$Patient_Time, colour="black", size=3, nudge_x=0.1, nudge_y=0.1, max.overlaps = Inf) +  
  xlab("PC1 (29%)")+
  ylab("PC2(10%)")+ 
  theme(axis.text=element_text(size=14, colour="black"),axis.title=element_text(size=14, face="bold", colour="black")) + 
  theme(panel.background = element_rect(fill = "white", colour = 'black'))+
  ggtitle("Study1_PCA_plot_without_outliers_Wagner_Study")
#dev.off()


##Differential gene expression analysis - Paired data
S1_subject_out <- study1_groups$Patient_Combined[match(colnames(study1_DGE_out), study1_groups$Raw.count.file)] #pull the patient and samples combined according to outlier removal adjustment
S1_subject_out <- factor(S1_subject_out)
S1_subject_out <- factor(S1_subject_out, levels = unique(S1_subject_out)) #get the subject data

## Time 
S1_Time_out <- study1_groups$Time[match(colnames(study1_DGE_out), study1_groups$Raw.count.file)]
S1_Time_out <- factor(S1_Time_out)
S1_Time_out <- factor(S1_Time_out, levels = unique(S1_Time_out))

S1_paired_data_out <- data.frame(Samplenames = colnames(study1_DGE_out$counts), Time = S1_Time_out, Subject = S1_subject_out)

## create a design matrix for the time comparison
S1_design_out <- model.matrix(~0+Time+Subject, data = S1_paired_data_out)
colnames(S1_design_out) <- gsub("Time", "", colnames(S1_design_out))
colnames(S1_design_out) <- gsub("Subject", "", colnames(S1_design_out))
rownames(S1_design_out) <- S1_paired_data_out$Samplenames

## contrast matrix for time comparison
S1_contr.mtrx_out <- makeContrasts(
  Time = post - pre,
  levels = S1_design_out)
S1_contr.mtrx_out

## run voom
S1_v_out <- voom(study1_DGE_out, S1_design_out, plot = T)
S1_fit_out <- lmFit(S1_v_out, S1_design_out)
S1_fit_out <- contrasts.fit(S1_fit_out, S1_contr.mtrx_out)
S1_efit_out <- eBayes(S1_fit_out)

PostvsPre_ALL_out <- topTable(S1_efit_out, n = Inf, adjust.method = "fdr", sort.by = "P") # no fdr threshold was given to this data - contains all genes
PostvsPre_FDR_0.05_out <- topTable(S1_efit_out, n = Inf, adjust.method = "fdr", sort.by = "P", p.value = 0.05) #found 5129 genes at FDR 0.05
PostvsPre_FDR_0.01_out <- topTable(S1_efit_out, n = Inf, adjust.method = "fdr", sort.by = "P", p.value = 0.01) #found 3442 genes at FDR 0.01

#wb <- loadWorkbook("/Users/JavvadP/Downloads/Wagner/out_limma_study1_paired.xlsx")
#addWorksheet(wb, "PostvsPre_ALL")
#writeData(wb, "PostvsPre_ALL",PostvsPre_ALL)
#saveWorkbook(wb, "/Users/JavvadP/Downloads/Wagner/out_limma_study1_paired.xlsx", overwrite = T)

## Differential gene expression analysis - Unpaired data
Group.Time_out <- factor(paste(study1_groups_out$Groups, S1_Time_out, sep = "_"))
Group.Time_out <- factor(Group.Time_out, levels = unique(Group.Time_out))
S1_group_out <- factor(study1_groups_out$Groups)
S1_group_out <- factor(S1_group_out, levels = unique(S1_group_out))

S1_unpaired_data_out <- data.frame(Samplenames = colnames(study1_DGE_out$counts), 
                        Time = S1_Time_out, Group = S1_group_out, Group.Time = Group.Time_out)


#design matrix for unpaired analysis
S1_design1_out <- model.matrix(~0+Group.Time_out, data = S1_unpaired_data_out)
colnames(S1_design1_out) <- gsub("Group.Time_out", "", colnames(S1_design1_out))
rownames(S1_design1_out) <- S1_unpaired_data_out$Samplenames

# contrast matrix for unpaired analysis
S1_contr.mtrx1_out <- makeContrasts(
  TNFvsTCZ = (TNF_pre + TNF_post)/2 - (TCZ_pre + TCZ_post)/2,
  TNFvsTCZ_pre = TNF_pre - TCZ_pre,
  TNFvsTCZ_post = TNF_post - TCZ_post,
  TNFvsTCZ_Intxn = (TNF_pre - TNF_post) - (TCZ_pre - TCZ_post),
  levels = S1_design1_out)

## Run voom
S1_v1_out <- voom(study1_DGE_out, S1_design1_out, plot = T)
S1_fit1_out <- lmFit(S1_v1_out, S1_design1_out)
S1_fit1_out <- contrasts.fit(S1_fit1_out, S1_contr.mtrx1_out)
S1_efit1_out <- eBayes(S1_fit1_out)

TNFvsTCZ_ALL_out <- topTable(S1_efit1_out, n = Inf, coef = 1, adjust.method = "fdr", sort.by = "P") # no fdr threshold was given to this data
TNFvsTCZ_pre_out <- topTable(S1_efit1_out, n = Inf, coef = 2, adjust.method = "fdr", sort.by = "P") # no fdr threshold was given to this data
TNFvsTCZ_post_out <- topTable(S1_efit1_out, n = Inf, coef = 3, adjust.method = "fdr", sort.by = "P") # no fdr threshold was given to this data
TNFvsTCZ_Intxn_out <- topTable(S1_efit1_out, n = Inf, coef = 4, adjust.method = "fdr", sort.by = "P") # no fdr threshold was given to this data

topTable(S1_efit1_out, n= Inf, coef = 1, adjust.method = "fdr", sort.by = "P", p.value = 0.05) # no genes were found at FDR 0.05 for TNFvsTCZ comparison
topTable(S1_efit1_out, n = Inf, coef = 2, adjust.method = "fdr", sort.by = "P", p.value = 0.05) # no genes were found at FDR 0.05 for TNFvsTCZ_pre comparison
topTable(S1_efit1_out, n = Inf, coef = 3, adjust.method = "fdr", sort.by = "P", p.value = 0.05) # no genes were found at FDR 0.05 for TNFvsTCZ_post comparison
topTable(S1_efit1_out, n = Inf, coef = 4, adjust.method = "fdr", sort.by = "P", p.value = 0.05) # no genes were found at FDR 0.05 for TNFvsTCZ_Intxn comparison

topTable(S1_efit1_out, n= Inf, coef = 1, adjust.method = "fdr", sort.by = "P", p.value = 0.1) # no genes were found at FDR 0.1 for TNFvsTCZ comparison
topTable(S1_efit1_out, n = Inf, coef = 2, adjust.method = "fdr", sort.by = "P", p.value = 0.1) # no genes were found at FDR 0.1 for TNFvsTCZ_pre comparison
topTable(S1_efit1_out, n = Inf, coef = 3, adjust.method = "fdr", sort.by = "P", p.value = 0.1) # no genes were found at FDR 0.1 for TNFvsTCZ_post comparison
topTable(S1_efit1_out, n = Inf, coef = 4, adjust.method = "fdr", sort.by = "P", p.value = 0.1) # no genes were found at FDR 0.1 for TNFvsTCZ_Intxn comparison

topTable(S1_efit1_out, n= Inf, coef = 1, adjust.method = "fdr", sort.by = "P", p.value = 0.2) # no genes were found at FDR 0.2 for TNFvsTCZ comparison
topTable(S1_efit1_out, n = Inf, coef = 2, adjust.method = "fdr", sort.by = "P", p.value = 0.2) # no genes were found at FDR 0.2 for TNFvsTCZ_pre comparison
topTable(S1_efit1_out, n = Inf, coef = 3, adjust.method = "fdr", sort.by = "P", p.value = 0.2) # no genes were found at FDR 0.2 for TNFvsTCZ_post comparison
topTable(S1_efit1_out, n = Inf, coef = 4, adjust.method = "fdr", sort.by = "P", p.value = 0.2) # no genes were found at FDR 0.2 for TNFvsTCZ_Intxn comparison



##--------------------------------------------------------Volcano Plots and Heatmaps---------------------------------

## Volcano plot for paired analysis
#jpeg("VolcanoPlot_PostvsPre_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot(PostvsPre_ALL_out, aes(x= logFC ,
                             y=-log10(adj.P.Val), 
                             color= -log10(adj.P.Val)))+
  geom_point(size=2)+
  geom_point(data=subset(PostvsPre_ALL_out %>% slice_min(adj.P.Val, n=20)),
             aes(x= logFC, 
                 y=-log10(adj.P.Val), 
                 color= -log10(adj.P.Val)))+
  scale_color_viridis(option="rocket", direction=1)+
  geom_text_repel(data=subset(PostvsPre_ALL_out %>% slice_min(adj.P.Val, n=20)),
                  aes(x=logFC, 
                      y=-log10(adj.P.Val), 
                      color= -log10(adj.P.Val),label = gene_name), 
                  color = "black", size=4, box.padding = unit(0.2, "lines"),segment.color="gray25",
                  segment.size=0.25, point.padding = unit(0.2, "lines"), max.overlaps = Inf)+
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(1,-1),linetype="dashed")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text (size = 10, face ="bold", colour = "black"))+
  theme(axis.title = element_text (size = 18, face = "bold", colour = "black")) +
  theme(strip.text = element_text(size=10, face="bold", color = "black"))+
  ggtitle("volcano_Plot_Top20genes_PostvsPre")+
  theme(plot.title=element_text(face="bold"))+
  xlab("log2FC (Post vs. Pre)") +
  ylab("-log10 FDR (Post vs. Pre)")+
  theme(legend.position="bottom") #to maximize graph space by removing legends
#dev.off()

## Volcano plot for unpaired analysis - TNFvsTCZ
#jpeg("volcanoplot_TNFvsTCZ_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot(TNFvsTCZ_ALL_out, aes(x= logFC ,
                          y=-log10(P.Value), 
                          color= -log10(P.Value)))+
  geom_point(size=2)+
  geom_point(data=subset(TNFvsTCZ_ALL_out %>% slice_min(P.Value, n=20)),
             aes(x= logFC, 
                 y=-log10(P.Value), 
                 color= -log10(P.Value)))+
  scale_color_viridis(option="rocket", direction=1)+
  geom_text_repel(data=subset(TNFvsTCZ_ALL_out %>% slice_min(P.Value, n=20)),
                  aes(x=logFC, 
                      y=-log10(P.Value), 
                      color= -log10(P.Value),label = gene_name), 
                  color = "black", size=4, box.padding = unit(0.2, "lines"),segment.color="gray25",
                  segment.size=0.25, point.padding = unit(0.2, "lines"), max.overlaps = Inf)+
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(1,-1),linetype="dashed")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text (size = 10, face ="bold", colour = "black"))+
  theme(axis.title = element_text (size = 18, face = "bold", colour = "black")) +
  theme(strip.text = element_text(size=10, face="bold", color = "black"))+
  ggtitle("volcano_Plot_Top20genes_TNFvsTCZ")+
  theme(plot.title=element_text(face="bold"))+
  xlab("log2FC (TNF vs.TCZ)") +
  ylab("-log10 P.Val (TNF vs. TCZ)")+
  theme(legend.position="bottom") #to maximize graph space by removing legends
#dev.off()

## Volcano plot for unpaired analysis - TNFvsTCZ_pre
#jpeg("volcanoplot_TNFvsTCZ_pre.jpeg", width = 2500, height = 1800, res = 300)
ggplot(TNFvsTCZ_pre_out, aes(x= logFC ,
                         y=-log10(P.Value), 
                         color= -log10(P.Value)))+
  geom_point(size=2)+
  geom_point(data=subset(TNFvsTCZ_pre_out %>% slice_min(P.Value, n=20)),
             aes(x= logFC, 
                 y=-log10(P.Value), 
                 color= -log10(P.Value)))+
  scale_color_viridis(option="rocket", direction=1)+
  geom_text_repel(data=subset(TNFvsTCZ_pre_out %>% slice_min(P.Value, n=20)),
                  aes(x=logFC, 
                      y=-log10(P.Value), 
                      color= -log10(P.Value),label = gene_name), 
                  color = "black", size=4, box.padding = unit(0.2, "lines"),segment.color="gray25",
                  segment.size=0.25, point.padding = unit(0.2, "lines"), max.overlaps = Inf)+
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(1,-1),linetype="dashed")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text (size = 10, face ="bold", colour = "black"))+
  theme(axis.title = element_text (size = 18, face = "bold", colour = "black")) +
  theme(strip.text = element_text(size=10, face="bold", color = "black"))+
  ggtitle("volcano_Plot_Top20genes_TNFvsTCZ_pre")+
  theme(plot.title=element_text(face="bold"))+
  xlab("log2FC (TNF_pre vs.TCZ_pre)") +
  ylab("-log10 P.Val (TNF_pre vs. TCZ_pre)")+
  theme(legend.position="bottom") #to maximize graph space by removing legends
#dev.off()


## Volcano plot for unpaired analysis - TNFvsTCZ_post
#jpeg("volcanoplot_TNFvsTCZ_post.jpeg", width = 2500, height = 1800, res = 300)
ggplot(TNFvsTCZ_post_out, aes(x= logFC ,
                         y=-log10(P.Value), 
                         color= -log10(P.Value)))+
  geom_point(size=2)+
  geom_point(data=subset(TNFvsTCZ_post_out %>% slice_min(P.Value, n=20)),
             aes(x= logFC, 
                 y=-log10(P.Value), 
                 color= -log10(P.Value)))+
  scale_color_viridis(option="rocket", direction=1)+
  geom_text_repel(data=subset(TNFvsTCZ_post_out %>% slice_min(P.Value, n=20)),
                  aes(x=logFC, 
                      y=-log10(P.Value), 
                      color= -log10(P.Value),label = gene_name), 
                  color = "black", size=4, box.padding = unit(0.2, "lines"),segment.color="gray25",
                  segment.size=0.25, point.padding = unit(0.2, "lines"), max.overlaps = Inf)+
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(1,-1),linetype="dashed")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text (size = 10, face ="bold", colour = "black"))+
  theme(axis.title = element_text (size = 18, face = "bold", colour = "black")) +
  theme(strip.text = element_text(size=10, face="bold", color = "black"))+
  ggtitle("volcano_Plot_Top20genes_TNFvsTCZ_post")+
  theme(plot.title=element_text(face="bold"))+
  xlab("log2FC (TNF_post vs.TCZ_post)") +
  ylab("-log10 P.Val (TNF_post vs. TCZ_post)")+
  theme(legend.position="bottom") #to maximize graph space by removing legends
#dev.off()

## Volcano plot for unpaired analysis - TNFvsTCZ_Intxn
#jpeg("volcanoplot_TNFvsTCZ_Intxn.jpeg", width = 3000, height = 2000, res = 300)
ggplot(TNFvsTCZ_Intxn_out, aes(x= logFC ,
                          y=-log10(P.Value), 
                          color= -log10(P.Value)))+
  geom_point(size=2)+
  geom_point(data=subset(TNFvsTCZ_Intxn_out %>% slice_min(P.Value, n=20)),
             aes(x= logFC, 
                 y=-log10(P.Value), 
                 color= -log10(P.Value)))+
  scale_color_viridis(option="rocket", direction=1)+
  geom_text_repel(data=subset(TNFvsTCZ_Intxn_out %>% slice_min(P.Value, n=20)),
                  aes(x=logFC, 
                      y=-log10(P.Value), 
                      color= -log10(P.Value),label = gene_name), 
                  color = "black", size=4, box.padding = unit(0.2, "lines"),segment.color="gray25",
                  segment.size=0.25, point.padding = unit(0.2, "lines"), max.overlaps = Inf)+
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(1,-1),linetype="dashed")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text (size = 10, face ="bold", colour = "black"))+
  theme(axis.title = element_text (size = 18, face = "bold", colour = "black")) +
  theme(strip.text = element_text(size=10, face="bold", color = "black"))+
  ggtitle("volcano_Plot_Top20genes_TNFvsTCZ_Intxn")+
  theme(plot.title=element_text(face="bold"))+
  xlab("log2FC (Interaction Effect)") +
  ylab("-log10 P.Val (Interaction_Effect)")+
  theme(legend.position="bottom") #to maximize graph space by removing legends
#dev.off()

## Heatmaps for paired analysis - PostvsPre
PostvsPre_ALL_out_filtered <- PostvsPre_ALL_out[PostvsPre_ALL_out$adj.P.Val < 0.05 & abs(PostvsPre_ALL_out$logFC) > 1, ]
i.1 <- S1_efit_out$genes$gene_name[which(S1_efit_out$genes$gene_name %in% PostvsPre_ALL_out_filtered$gene_name)] 

data_PostvsPre_ALL_out <- study1_logCPM_out_file[study1_logCPM_out_file$Gene %in% i.1, -1]
rownames(data_PostvsPre_ALL_out) <- i.1

samples.time <- paste(colnames(data_PostvsPre_ALL_out), S1_Time_out, sep = "_")
colnames(data_PostvsPre_ALL_out) <- samples.time

pre_samples <- grep("pre", colnames(data_PostvsPre_ALL_out), value = TRUE)  # Extract pre samples
post_samples <- grep("post", colnames(data_PostvsPre_ALL_out), value = T) #extract post samples

# Reorder the columns such that pre and post samples are grouped
S1_ordered_samples <- c(pre_samples, post_samples)

#color palette
palette <- colorRampPalette(c("blue", "white", "red"))(100)

data_PostvsPre_ALL_out <- data_PostvsPre_ALL_out[, S1_ordered_samples]

# create an annotation dataframe for the samples
annot_col_PrevsPost <- data.frame(Time = substr(colnames(data_PostvsPre_ALL_out),6,9))
rownames(annot_col_PrevsPost) <- colnames(data_PostvsPre_ALL_out)


annot_colors_PrevsPost <- list(
  Time = c(pre = "#1B9E77", post = "#D95F02")
)

data_PostvsPre_ALL_out <- data_PostvsPre_ALL_out %>% slice(1:100)

#jpeg("Heatmap_PostvsPre.jpeg", width = 3000, height = 3500, res = 300)
#generate heatmap
pheatmap(data_PostvsPre_ALL_out, show_rownames = T, scale = "row", cluster_cols= F, cluster_rows = T, color = palette, border_color = NA,
         fontsize_row = 8 , fontsize_col = 10,annotation_col = annot_col_PrevsPost, annotation_colors = annot_colors_PrevsPost, 
         annotation_legend = T, main = "Topgenes_PrevsPost_FDR_0.05")
#dev.off()


## Heatmap for unpaired analysis - TNFvsTCZ
TNFvsTCZ_ALL_out_filtered <- TNFvsTCZ_ALL_out[TNFvsTCZ_ALL_out$P.Value < 0.01 & abs(TNFvsTCZ_ALL_out$logFC) > 1, ]
i.2 <- S1_efit1_out$genes$gene_name[which(S1_efit1_out$genes$gene_name %in% TNFvsTCZ_ALL_out_filtered$gene_name)]

data_TNFvsTCZ_ALL_out <- study1_logCPM_out_file[study1_logCPM_out_file$Gene %in% i.2, -1]
rownames(data_TNFvsTCZ_ALL_out) <- i.2

colnames(data_TNFvsTCZ_ALL_out) <- samples.time

data_TNFvsTCZ_ALL_out <- data_TNFvsTCZ_ALL_out[, S1_ordered_samples]

# create an annotation dataframe for the samples
annot_col_TNFvsTCZ <- data.frame(Group = Group.Time_out)
rownames(annot_col_TNFvsTCZ) <- samples.time

# re arrange the rownames so that pre samples are followed by post samples and the group column values are adjusted accordingly
annot_col_TNFvsTCZ <- annot_col_TNFvsTCZ[order(grepl("post", rownames(annot_col_TNFvsTCZ)), rownames(annot_col_TNFvsTCZ)), , drop = FALSE]

annot_col_TNFvsTCZ$Group <- sub("_pre|_post", "", annot_col_TNFvsTCZ$Group) # remove pre and post suffix from the group column values


annot_colors_TNFvsTCZ <- list(
  Group = c(TNF = "#1B9E77", TCZ = "#D95F02")
)

#jpeg("Heatmap_TNFvsTCZ.jpeg", width = 2500, height = 1800, res = 300)
#generate heatmap
pheatmap(data_TNFvsTCZ_ALL_out, show_rownames = T, scale = "row", cluster_cols= F, cluster_rows = T, color = palette, border_color = NA,
         fontsize_row = 8 , fontsize_col = 10,annotation_col = annot_col_TNFvsTCZ, annotation_colors = annot_colors_TNFvsTCZ, 
         annotation_legend = T, main = "Topgenes_TNFvsTCZ_P.Val_0.01")
#dev.off()


## Heatmap for unpaired analysis - TNFvsTCZ_pre
TNFvsTCZ_pre_out_filtered <- TNFvsTCZ_pre_out[TNFvsTCZ_pre_out$P.Value < 0.01 & abs(TNFvsTCZ_pre_out$logFC) > 1, ] # filter the TNFvsTCZ_pre values based on the P.val < 0.01 and logfc > 1
i.3 <- S1_efit1_out$genes$gene_name[which(S1_efit1_out$genes$gene_name %in% TNFvsTCZ_pre_out_filtered$gene_name)] 

data_TNFvsTCZ_pre_out <- study1_logCPM_out_file[study1_logCPM_out_file$Gene %in% i.3, -1]
rownames(data_TNFvsTCZ_pre_out) <- i.3

colnames(data_TNFvsTCZ_pre_out) <- samples.time

data_TNFvsTCZ_pre_out <-  data_TNFvsTCZ_pre_out[, grepl("pre", colnames(data_TNFvsTCZ_pre_out))] # keep only the pre samples

# create an annotation dataframe for the samples
annot_col_TNFvsTCZ_pre <- data.frame(Group = Group.Time_out)
rownames(annot_col_TNFvsTCZ_pre) <- samples.time

# re arrange the rownames so that pre samples are followed by post samples and the group column values are adjusted accordingly
annot_col_TNFvsTCZ_pre <- annot_col_TNFvsTCZ_pre[order(grepl("post", rownames(annot_col_TNFvsTCZ_pre)), rownames(annot_col_TNFvsTCZ_pre)), , drop = FALSE]

annot_col_TNFvsTCZ_pre <- annot_col_TNFvsTCZ_pre[!grepl("post", rownames(annot_col_TNFvsTCZ_pre)), , drop = F] #remove the post samples related rows

annot_col_TNFvsTCZ_pre$Group <- sub("_pre", "", annot_col_TNFvsTCZ_pre$Group) # remove pre suffix from the group column values

annot_colors_TNFvsTCZ_pre <- list(
  Group = c(TNF = "#1B9E77", TCZ = "#D95F02")
)

#jpeg("Heatmap_TNFvsTCZ_pre.jpeg", width = 2500, height = 1800, res = 300)
#generate heatmap
pheatmap(data_TNFvsTCZ_pre_out, show_rownames = T, scale = "row", cluster_cols= F, cluster_rows = T, color = palette, border_color = NA,
         fontsize_row = 8 , fontsize_col = 10,annotation_col = annot_col_TNFvsTCZ_pre, annotation_colors = annot_colors_TNFvsTCZ_pre, 
         annotation_legend = T, main = "Topgenes_TNFvsTCZ_pre_P.Val_0.01")
#dev.off()

## Heatmap for unpaired analysis - TNFvsTCZ_post
TNFvsTCZ_post_out_filtered <- TNFvsTCZ_post_out[TNFvsTCZ_post_out$P.Value < 0.01 & abs(TNFvsTCZ_post_out$logFC) > 1, ] # filter the TNFvsTCZ_post values based on the P.val < 0.01 and logfc > 1
i.4 <- S1_efit1_out$genes$gene_name[which(S1_efit1_out$genes$gene_name %in% TNFvsTCZ_post_out_filtered$gene_name)] 

data_TNFvsTCZ_post_out <- study1_logCPM_out_file[study1_logCPM_out_file$Gene %in% i.4, -1]
rownames(data_TNFvsTCZ_post_out) <- i.4

colnames(data_TNFvsTCZ_post_out) <- samples.time

data_TNFvsTCZ_post_out <-  data_TNFvsTCZ_post_out[, grepl("post", colnames(data_TNFvsTCZ_post_out))] # keep only the post samples 

# create an annotation dataframe for the samples
annot_col_TNFvsTCZ_post <- data.frame(Group = Group.Time_out)
rownames(annot_col_TNFvsTCZ_post) <- samples.time

# re arrange the rownames so that pre samples are followed by post samples and the group column values are adjusted accordingly
annot_col_TNFvsTCZ_post <- annot_col_TNFvsTCZ_post[order(grepl("post", rownames(annot_col_TNFvsTCZ_post)), rownames(annot_col_TNFvsTCZ_post)), , drop = FALSE]

annot_col_TNFvsTCZ_post <- annot_col_TNFvsTCZ_post[!grepl("pre", rownames(annot_col_TNFvsTCZ_post)), , drop = F] #remove the pre samples related rows


annot_col_TNFvsTCZ_post$Group <- sub("_post", "", annot_col_TNFvsTCZ_post$Group) # remove post suffix from the group column values

annot_colors_TNFvsTCZ_post <- list(
  Group = c(TNF = "#1B9E77", TCZ = "#D95F02")
)

#jpeg("Heatmap_TNFvsTCZ_post.jpeg", width = 2500, height = 1800, res = 300)
#generate heatmap
pheatmap(data_TNFvsTCZ_post_out, show_rownames = T, scale = "row", cluster_cols= F, cluster_rows = T, color = palette, border_color = NA,
         fontsize_row = 8 , fontsize_col = 10,annotation_col = annot_col_TNFvsTCZ_post, annotation_colors = annot_colors_TNFvsTCZ_post, 
         annotation_legend = T, main = "Topgenes_TNFvsTCZ_post_P.Val_0.01")
#dev.off()


## Heatmap for unpaired analysis - Interaction
TNFvsTCZ_Intxn_out_filtered <- TNFvsTCZ_Intxn_out[TNFvsTCZ_Intxn_out$P.Value < 0.01 & abs(TNFvsTCZ_Intxn_out$logFC) > 1, ] # filter the TNFvsTCZ_post values based on the P.val < 0.01 and logfc > 1
i.5 <- S1_efit1_out$genes$gene_name[which(S1_efit1_out$genes$gene_name %in% TNFvsTCZ_Intxn_out_filtered$gene_name)] 

data_TNFvsTCZ_Intxn_out <- study1_logCPM_out_file[study1_logCPM_out_file$Gene %in% i.5, -1]
rownames(data_TNFvsTCZ_Intxn_out) <- i.5

colnames(data_TNFvsTCZ_Intxn_out) <- samples.time

data_TNFvsTCZ_Intxn_out <- data_TNFvsTCZ_Intxn_out[, S1_ordered_samples]# order the samples such that pre samples are followed by post samples

# create an annotation dataframe for the samples
annot_col_TNFvsTCZ_Intxn <- data.frame(Group = Group.Time_out)
rownames(annot_col_TNFvsTCZ_Intxn) <- samples.time

# re arrange the rownames so that pre samples are followed by post samples and the group column values are adjusted accordingly
annot_col_TNFvsTCZ_Intxn <- annot_col_TNFvsTCZ_Intxn[order(grepl("post", rownames(annot_col_TNFvsTCZ_Intxn)), rownames(annot_col_TNFvsTCZ_Intxn)), , drop = FALSE]

annot_col_TNFvsTCZ_Intxn$Time <- ifelse(               # extract timepoints from group column
  grepl("pre", annot_col_TNFvsTCZ_Intxn$Group), "pre",
  ifelse(grepl("post", annot_col_TNFvsTCZ_Intxn$Group), "post", NA)
)

annot_col_TNFvsTCZ_Intxn$Group <- sub("_pre|_post", "", annot_col_TNFvsTCZ_Intxn$Group) # remove pre and post suffix from the group column values

annot_colors_TNFvsTCZ_Intxn <- list(
  Group = c(TNF = "#1B9E77", TCZ = "#D95F02"),
  Time = c(pre = "darkblue", post = "darkred")
)

#jpeg("Heatmap_TNFvsTCZ_Intxn.jpeg", width = 2500, height = 1800, res = 300)
#generate heatmap
pheatmap(data_TNFvsTCZ_Intxn_out, show_rownames = T, scale = "row", cluster_cols= F, cluster_rows = T, color = palette, border_color = NA,
         fontsize_row = 8 , fontsize_col = 10,annotation_col = annot_col_TNFvsTCZ_Intxn, annotation_colors = annot_colors_TNFvsTCZ_Intxn, 
         annotation_legend = T, main = "Topgenes_TNFvsTCZ_Intxn_P.Val_0.01")
#dev.off()

##--------------------------------------Gene Set Enrichment Analysis - PostvsPre_ALL---------------------------
Hallmark <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt") # from MSigDB
Wiki <- gmtPathways("c2.cp.wikipathways.v2024.1.Hs.symbols.gmt") #from MSigDB
KEGG <- gmtPathways("c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt") #from MSigDB

# prepare ranked gene list for Post vs Pre
ranked_genes_PostvsPre <- PostvsPre_ALL_out$logFC
names(ranked_genes_PostvsPre) <- PostvsPre_ALL_out$gene_name

# Sort genes in decreasing order of logFC
ranked_genes_PostvsPre <- sort(ranked_genes_PostvsPre, decreasing = TRUE)

# View the ranked genes
head(ranked_genes_PostvsPre)

# Run fgsea with HPostvsPremark pathways 
FGSEA_hl_PostvsPre <- fgsea(pathways = Hallmark , # List of gene sets to check
                      stats = ranked_genes_PostvsPre,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_hl_PostvsPre$leadingEdge <- sapply(FGSEA_hl_PostvsPre$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_hl_PostvsPre, "Fgsea_Hallmark_PostvsPre.xlsx")

#select top pathways and create a table
topPathwaysUp_hl_PostvsPre<- FGSEA_hl_PostvsPre[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 hallmark pathways upregulated at pval 0.01
topPathwaysDown_hl_PostvsPre <- FGSEA_hl_PostvsPre[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 hallmark pathways downregulated at pval 0.01
topPathways_hl_PostvsPre <- c(topPathwaysUp_hl_PostvsPre, rev(topPathwaysDown_hl_PostvsPre))
plotGseaTable(Hallmark[topPathways_hl_PostvsPre], stats = ranked_genes_PostvsPre, fgseaRes = FGSEA_hl_PostvsPre, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_hl_PostvsPre <- collapsePathways(FGSEA_hl_PostvsPre[order(pval)][padj < 0.01], Hallmark, ranked_genes_PostvsPre, gseaParam = 1)

# plot the most significantly enriched pathway
plotEnrichment(Hallmark[[head(FGSEA_hl_PostvsPre[order(pval), ], 1)$pathway]],
               ranked_genes_PostvsPre) + 
  labs(title = head(FGSEA_hl_PostvsPre[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_hl_PostvsPre = Hallmark[topPathways_hl_PostvsPre]

# generate data for input into pathway images
df_hl_PostvsPre = gseaCurve(ranked_genes_PostvsPre, setlist_hl_PostvsPre, FGSEA_hl_PostvsPre)

#jpeg("HallmarkPathway_PostvsPre_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_hl_PostvsPre, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "HallmarkPathway_PostvsPre_ALL(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()


## Wiki pathways
# Run fgsea with Wiki pathways 
FGSEA_wk_PostvsPre <- fgsea(pathways = Wiki , # List of gene sets to check
                      stats = ranked_genes_PostvsPre,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_wk_PostvsPre$leadingEdge <- sapply(FGSEA_wk_PostvsPre$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_wk_PostvsPre, "Fgsea_Wiki_PostvsPre.xlsx")

#select top pathways and create a table
topPathwaysUp_wk_PostvsPre <- FGSEA_wk_PostvsPre[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 wiki pathways upregulated at pval 0.01
topPathwaysDown_wk_PostvsPre <- FGSEA_wk_PostvsPre[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 wiki pathways downregulated at pval 0.01
topPathways_wk_PostvsPre <- c(topPathwaysUp_wk_PostvsPre, rev(topPathwaysDown_wk_PostvsPre))
plotGseaTable(Wiki[topPathways_wk_PostvsPre], stats = ranked_genes_PostvsPre, fgseaRes = FGSEA_wk_PostvsPre, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_wk_PostvsPre <- collapsePathways(FGSEA_wk_PostvsPre[order(pval)][padj < 0.01], Wiki, ranked_genes_PostvsPre, gseaParam = 1)
mainPathways_wk_PostvsPre <- FGSEA_wk_PostvsPre[pathway %in% collapsedPathways_wk_PostvsPre$mainPathways][order(-NES), pathway]
plotGseaTable(Wiki[mainPathways_wk_PostvsPre], ranked_genes_PostvsPre, FGSEA_wk_PostvsPre, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(Wiki[[head(FGSEA_wk_PostvsPre[order(pval), ], 1)$pathway]],
               ranked_genes_PostvsPre) + 
  labs(title = head(FGSEA_wk_PostvsPre[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_wk_PostvsPre = Wiki[topPathways_wk_PostvsPre]

# generate data for input into pathway images
df_wk_PostvsPre = gseaCurve(ranked_genes_PostvsPre, setlist_wk_PostvsPre, FGSEA_wk_PostvsPre)

#jpeg("WikiPathway_PostvsPre_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_wk_PostvsPre, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "WikiPathway_PostvsPre_ALL(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

## KEGG pathway
# Run fgsea with KEGG pathways 
FGSEA_kg_PostvsPre <- fgsea(pathways = KEGG , # List of gene sets to check
                      stats = ranked_genes_PostvsPre,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_kg_PostvsPre$leadingEdge <- sapply(FGSEA_kg_PostvsPre$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_kg_PostvsPre, "Fgsea_KEGG_PostvsPre.xlsx")

#select top pathways and create a table
topPathwaysUp_kg_PostvsPre <- FGSEA_kg_PostvsPre[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 kegg pathways upregulated at pval 0.01
topPathwaysDown_kg_PostvsPre <- FGSEA_kg_PostvsPre[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 kegg pathways downregulated at pval 0.01
topPathways_kg_PostvsPre <- c(topPathwaysUp_kg_PostvsPre, rev(topPathwaysDown_kg_PostvsPre))
plotGseaTable(KEGG[topPathways_kg_PostvsPre], stats = ranked_genes_PostvsPre, fgseaRes = FGSEA_kg_PostvsPre, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_kg_PostvsPre <- collapsePathways(FGSEA_kg_PostvsPre[order(pval)][padj < 0.01], KEGG, ranked_genes_PostvsPre, gseaParam = 1)
mainPathways_kg_PostvsPre <- FGSEA_kg_PostvsPre[pathway %in% collapsedPathways_kg_PostvsPre$mainPathways][order(-NES), pathway]
plotGseaTable(KEGG[mainPathways_kg_PostvsPre], ranked_genes_PostvsPre, FGSEA_kg_PostvsPre, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(KEGG[[head(FGSEA_kg_PostvsPre[order(pval), ], 1)$pathway]],
               ranked_genes_PostvsPre) + 
  labs(title = head(FGSEA_kg_PostvsPre[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_kg_PostvsPre = KEGG[topPathways_kg_PostvsPre]

# generate data for input into pathway images
df_kg_PostvsPre = gseaCurve(ranked_genes_PostvsPre, setlist_kg_PostvsPre, FGSEA_kg_PostvsPre)

#jpeg("KEGGPathway_PostvsPre_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_kg_PostvsPre, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "KEGGPathway_PostvsPre_ALL(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

#FGSEA_hl_PostvsPre <- FGSEA_hl_PostvsPre %>% mutate(Database = "Hallmark")
#FGSEA_wk_PostvsPre <- FGSEA_wk_PostvsPre %>% mutate(Database = "Wiki")
#FGSEA_kg_PostvsPre <- FGSEA_kg_PostvsPre %>% mutate(Database = "KEGG")

#combine_PostvsPre <- bind_rows(FGSEA_hl_PostvsPre, FGSEA_wk_PostvsPre, FGSEA_kg_PostvsPre)
#combine_PostvsPre <- combine_PostvsPre %>% select(Database, everything())

#wb <- loadWorkbook("out_fgsea_study1.xlsx")
#addWorksheet(wb, "PostvsPre_ALL")
#writeData(wb, "PostvsPre_ALL", combine_PostvsPre)
#saveWorkbook(wb, "out_fgsea_Study1.xlsx", overwrite = T)

##--------------------------------------Barplots- FGSEA - PostvsPre_ALL------------------------------------------
top_hl_PostvsPre <- FGSEA_hl_PostvsPre[pathway %in% topPathways_hl_PostvsPre][order(pval)] #pull the top pathways from hallmark db

#Barplot for hallmark pathway - PostvsPre

#jpeg("FGSEA_Hallmark_PostvsPre_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_hl_PostvsPre) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Hallmark Pathway - PostvsPre_ALL")
#dev.off()


## Wiki pathways
top_wk_PostvsPre <- FGSEA_wk_PostvsPre[pathway %in% topPathways_wk_PostvsPre][order(pval)] #pull the top pathways from wiki db

#Barplot for wiki pathway - PostvsPre

#jpeg("FGSEA_Wiki_PostvsPre_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_wk_PostvsPre) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Wiki Pathway - PostvsPre_ALL")
#dev.off()

## KEGG pathways
top_kg_PostvsPre <- FGSEA_kg_PostvsPre[pathway %in% topPathways_kg_PostvsPre][order(pval)] #pull the top pathways from KEGG db

#Barplot for KEGG pathway - PostvsPre

#jpeg("FGSEA_KEGG_PostvsPre_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_kg_PostvsPre) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="KEGG Pathway - PostvsPre_ALL")
#dev.off()

##--------------------------------------Gene Set Enrichment Analysis - TNFvsTCZ_ALL---------------------------
Hallmark <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt") # from MSigDB
Wiki <- gmtPathways("c2.cp.wikipathways.v2024.1.Hs.symbols.gmt") #from MSigDB
KEGG <- gmtPathways("c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt") #from MSigDB

# prepare ranked gene list for Pre vs Post
ranked_genes_ALL <- TNFvsTCZ_ALL_out$logFC
names(ranked_genes_ALL) <- TNFvsTCZ_ALL_out$gene_name

# Sort genes in decreasing order of logFC
ranked_genes_ALL <- sort(ranked_genes_ALL, decreasing = TRUE)

# View the ranked genes
head(ranked_genes_ALL)

# Run fgsea with Hallmark pathways 
FGSEA_hl_ALL <- fgsea(pathways = Hallmark , # List of gene sets to check
                      stats = ranked_genes_ALL,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_hl_ALL$leadingEdge <- sapply(FGSEA_hl_ALL$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_hl_ALL, "Fgsea_Hallmark_ALL.xlsx")

#select top pathways and create a table
topPathwaysUp_hl_ALL <- FGSEA_hl_ALL[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 hallmark pathways upregulated at pval 0.01
topPathwaysDown_hl_ALL <- FGSEA_hl_ALL[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 hallmark pathways downregulated at pval 0.01
topPathways_hl_ALL <- c(topPathwaysUp_hl_ALL, rev(topPathwaysDown_hl_ALL))
plotGseaTable(Hallmark[topPathways_hl_ALL], stats = ranked_genes_ALL, fgseaRes = FGSEA_hl_ALL, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_hl_ALL <- collapsePathways(FGSEA_hl_ALL[order(pval)][padj < 0.01], Hallmark, ranked_genes_ALL, gseaParam = 1)

# plot the most significantly enriched pathway
plotEnrichment(Hallmark[[head(FGSEA_hl_ALL[order(pval), ], 1)$pathway]],
               ranked_genes_ALL) + 
  labs(title = head(FGSEA_hl_ALL[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_hl_ALL = Hallmark[topPathways_hl_ALL]

# generate data for input into pathway images
df_hl_ALL = gseaCurve(ranked_genes_ALL, setlist_hl_ALL, FGSEA_hl_ALL)

#jpeg("HallmarkPathway_TNFvsTCZ_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_hl_ALL, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "HallmarkPathway_TNFvsTCZ_ALL(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()


## Wiki pathways
# Run fgsea with Wiki pathways 
FGSEA_wk_ALL <- fgsea(pathways = Wiki , # List of gene sets to check
                      stats = ranked_genes_ALL,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_wk_ALL$leadingEdge <- sapply(FGSEA_wk_ALL$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_wk_ALL, "Fgsea_Wiki_ALL.xlsx")

#select top pathways and create a table
topPathwaysUp_wk_ALL <- FGSEA_wk_ALL[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 wiki pathways upregulated at pval 0.01
topPathwaysDown_wk_ALL <- FGSEA_wk_ALL[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 wiki pathways downregulated at pval 0.01
topPathways_wk_ALL <- c(topPathwaysUp_wk_ALL, rev(topPathwaysDown_wk_ALL))
plotGseaTable(Wiki[topPathways_wk_ALL], stats = ranked_genes_ALL, fgseaRes = FGSEA_wk_ALL, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_wk_ALL <- collapsePathways(FGSEA_wk_ALL[order(pval)][padj < 0.01], Wiki, ranked_genes_ALL, gseaParam = 1)
mainPathways_wk_ALL <- FGSEA_wk_ALL[pathway %in% collapsedPathways_wk_ALL$mainPathways][order(-NES), pathway]
plotGseaTable(Wiki[mainPathways_wk_ALL], ranked_genes_ALL, FGSEA_wk_ALL, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(Wiki[[head(FGSEA_wk_ALL[order(pval), ], 1)$pathway]],
               ranked_genes_ALL) + 
  labs(title = head(FGSEA_wk_ALL[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_wk_ALL = Wiki[topPathways_wk_ALL]

# generate data for input into pathway images
df_wk_ALL = gseaCurve(ranked_genes_ALL, setlist_wk_ALL, FGSEA_wk_ALL)

#jpeg("WikiPathway_TNFvsTCZ_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_wk_ALL, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "WikiPathway_TNFvsTCZ_ALL(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

## KEGG pathway
# Run fgsea with KEGG pathways 
FGSEA_kg_ALL <- fgsea(pathways = KEGG , # List of gene sets to check
                      stats = ranked_genes_ALL,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_kg_ALL$leadingEdge <- sapply(FGSEA_kg_ALL$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_kg_ALL, "Fgsea_KEGG_ALL.xlsx")

#select top pathways and create a table
topPathwaysUp_kg_ALL <- FGSEA_kg_ALL[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 kegg pathways upregulated at pval 0.01
topPathwaysDown_kg_ALL <- FGSEA_kg_ALL[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 kegg pathways downregulated at pval 0.01
topPathways_kg_ALL <- c(topPathwaysUp_kg_ALL, rev(topPathwaysDown_kg_ALL))
plotGseaTable(KEGG[topPathways_kg_ALL], stats = ranked_genes_ALL, fgseaRes = FGSEA_kg_ALL, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_kg_ALL <- collapsePathways(FGSEA_kg_ALL[order(pval)][padj < 0.01], KEGG, ranked_genes_ALL, gseaParam = 1)
mainPathways_kg_ALL <- FGSEA_kg_ALL[pathway %in% collapsedPathways_kg_ALL$mainPathways][order(-NES), pathway]
plotGseaTable(KEGG[mainPathways_kg_ALL], ranked_genes_ALL, FGSEA_kg_ALL, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(KEGG[[head(FGSEA_kg_ALL[order(pval), ], 1)$pathway]],
               ranked_genes_ALL) + 
  labs(title = head(FGSEA_kg_ALL[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_kg_ALL = KEGG[topPathways_kg_ALL]

# generate data for input into pathway images
df_kg_ALL = gseaCurve(ranked_genes_ALL, setlist_kg_ALL, FGSEA_kg_ALL)

#jpeg("KEGGPathway_TNFvsTCZ_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_kg_ALL, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "KEGGPathway_TNFvsTCZ_ALL(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

#FGSEA_hl_ALL <- FGSEA_hl_ALL %>% mutate(Database = "Hallmark")
#FGSEA_wk_ALL <- FGSEA_wk_ALL %>% mutate(Database = "Wiki")
#FGSEA_kg_ALL <- FGSEA_kg_ALL %>% mutate(Database = "KEGG")

#combine_ALL <- bind_rows(FGSEA_hl_ALL, FGSEA_wk_ALL, FGSEA_kg_ALL)
#combine_ALL <- combine_ALL %>% select(Database, everything())

#wb <- loadWorkbook("out_fgsea_study1.xlsx")
#addWorksheet(wb, "TNFvsTCZ_ALL")
#writeData(wb, "TNFvsTCZ_ALL", combine_ALL)
#saveWorkbook(wb, "out_fgsea_Study1.xlsx", overwrite = T)


##--------------------------------------Barplots- FGSEA - TNFvsTCZ_ALL------------------------------------------
top_hl_ALL <- FGSEA_hl_ALL[pathway %in% topPathways_hl_ALL][order(pval)] #pull the top pathways from hallmark db

#Barplot for hallmark pathway - TNFvsTCZ_ALL

#jpeg("FGSEA_Hallmark_TNFvsTCZ_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_hl_ALL) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Hallmark Pathway - TNFvsTCZ_ALL")
#dev.off()


## Wiki pathways
top_wk_ALL <- FGSEA_wk_ALL[pathway %in% topPathways_wk_ALL][order(pval)] #pull the top pathways from wiki db

#Barplot for wiki pathway - TNFvsTCZ_ALL

#jpeg("FGSEA_Wiki_TNFvsTCZ_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_wk_ALL) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Wiki Pathway - TNFvsTCZ_ALL")
#dev.off()

## KEGG pathways
top_kg_ALL <- FGSEA_kg_ALL[pathway %in% topPathways_kg_ALL][order(pval)] #pull the top pathways from KEGG db

#Barplot for KEGG pathway - TNFvsTCZ_ALL

#jpeg("FGSEA_KEGG_TNFvsTCZ_ALL.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_kg_ALL) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="KEGG Pathway - TNFvsTCZ_ALL")
#dev.off()



##--------------------------------------Gene Set Enrichment Analysis - TNFvsTCZ_Pre---------------------------
Hallmark <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt") # from MSigDB
Wiki <- gmtPathways("c2.cp.wikipathways.v2024.1.Hs.symbols.gmt") #from MSigDB
KEGG <- gmtPathways("c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt") #from MSigDB

# prepare ranked gene list for Pre vs Post
ranked_genes_Pre <- TNFvsTCZ_pre_out$logFC
names(ranked_genes_Pre) <- TNFvsTCZ_pre_out$gene_name

# Sort genes in decreasing order of logFC
ranked_genes_Pre <- sort(ranked_genes_Pre, decreasing = TRUE)

# View the ranked genes
head(ranked_genes_Pre)

# Run fgsea with Hallmark pathways 
FGSEA_hl_Pre <- fgsea(pathways = Hallmark , # List of gene sets to check
                      stats = ranked_genes_Pre,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_hl_Pre$leadingEdge <- sapply(FGSEA_hl_Pre$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_hl_Pre, "Fgsea_Hallmark_Pre.xlsx")

#select top pathways and create a table
topPathwaysUp_hl_Pre <- FGSEA_hl_Pre[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 hallmark pathways upregulated at pval 0.01
topPathwaysDown_hl_Pre <- FGSEA_hl_Pre[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 hallmark pathways downregulated at pval 0.01
topPathways_hl_Pre <- c(topPathwaysUp_hl_Pre, rev(topPathwaysDown_hl_Pre))
plotGseaTable(Hallmark[topPathways_hl_Pre], stats = ranked_genes_Pre, fgseaRes = FGSEA_hl_Pre, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_hl_Pre <- collapsePathways(FGSEA_hl_Pre[order(pval)][padj < 0.01], Hallmark, ranked_genes_Pre, gseaParam = 1)

# plot the most significantly enriched pathway
plotEnrichment(Hallmark[[head(FGSEA_hl_Pre[order(pval), ], 1)$pathway]],
               ranked_genes_Pre) + 
  labs(title = head(FGSEA_hl_Pre[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_hl_Pre = Hallmark[topPathways_hl_Pre]

# generate data for input into pathway images
df_hl_Pre = gseaCurve(ranked_genes_Pre, setlist_hl_Pre, FGSEA_hl_Pre)

jpeg("HallmarkPathway_TNFvsTCZ_Pre.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_hl_Pre, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "HallmarkPathway_TNFvsTCZ_Pre(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
dev.off()


## Wiki pathways
# Run fgsea with Wiki pathways 
FGSEA_wk_Pre <- fgsea(pathways = Wiki , # List of gene sets to check
                      stats = ranked_genes_Pre,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_wk_Pre$leadingEdge <- sapply(FGSEA_wk_Pre$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_wk_Pre, "Fgsea_Wiki_Pre.xlsx")

#select top pathways and create a table
topPathwaysUp_wk_Pre <- FGSEA_wk_Pre[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 wiki pathways upregulated at pval 0.01
topPathwaysDown_wk_Pre <- FGSEA_wk_Pre[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 wiki pathways downregulated at pval 0.01
topPathways_wk_Pre <- c(topPathwaysUp_wk_Pre, rev(topPathwaysDown_wk_Pre))
plotGseaTable(Wiki[topPathways_wk_Pre], stats = ranked_genes_Pre, fgseaRes = FGSEA_wk_Pre, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_wk_Pre <- collapsePathways(FGSEA_wk_Pre[order(pval)][padj < 0.01], Wiki, ranked_genes_Pre, gseaParam = 1)
mainPathways_wk_Pre <- FGSEA_wk_Pre[pathway %in% collapsedPathways_wk_Pre$mainPathways][order(-NES), pathway]
plotGseaTable(Wiki[mainPathways_wk_Pre], ranked_genes_Pre, FGSEA_wk_Pre, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(Wiki[[head(FGSEA_wk_Pre[order(pval), ], 1)$pathway]],
               ranked_genes_Pre) + 
  labs(title = head(FGSEA_wk_Pre[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_wk_Pre = Wiki[topPathways_wk_Pre]

# generate data for input into pathway images
df_wk_Pre = gseaCurve(ranked_genes_Pre, setlist_wk_Pre, FGSEA_wk_Pre)

#jpeg("WikiPathway_TNFvsTCZ_Pre.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_wk_Pre, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "WikiPathway_TNFvsTCZ_Pre(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

## KEGG pathway
# Run fgsea with KEGG pathways 
FGSEA_kg_Pre <- fgsea(pathways = KEGG , # List of gene sets to check
                      stats = ranked_genes_Pre,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_kg_Pre$leadingEdge <- sapply(FGSEA_kg_Pre$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_kg_Pre, "Fgsea_KEGG_Pre.xlsx")

#select top pathways and create a table
topPathwaysUp_kg_Pre <- FGSEA_kg_Pre[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 kegg pathways upregulated at pval 0.01
topPathwaysDown_kg_Pre <- FGSEA_kg_Pre[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 kegg pathways downregulated at pval 0.01
topPathways_kg_Pre <- c(topPathwaysUp_kg_Pre, rev(topPathwaysDown_kg_Pre))
plotGseaTable(KEGG[topPathways_kg_Pre], stats = ranked_genes_Pre, fgseaRes = FGSEA_kg_Pre, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_kg_Pre <- collapsePathways(FGSEA_kg_Pre[order(pval)][padj < 0.01], KEGG, ranked_genes_Pre, gseaParam = 1)
mainPathways_kg_Pre <- FGSEA_kg_Pre[pathway %in% collapsedPathways_kg_Pre$mainPathways][order(-NES), pathway]
plotGseaTable(KEGG[mainPathways_kg_Pre], ranked_genes_Pre, FGSEA_kg_Pre, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(KEGG[[head(FGSEA_kg_Pre[order(pval), ], 1)$pathway]],
               ranked_genes_Pre) + 
  labs(title = head(FGSEA_kg_Pre[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_kg_Pre = KEGG[topPathways_kg_Pre]

# generate data for input into pathway images
df_kg_Pre = gseaCurve(ranked_genes_Pre, setlist_kg_Pre, FGSEA_kg_Pre)

#jpeg("KEGGPathway_TNFvsTCZ_Pre.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_kg_Pre, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "KEGGPathway_TNFvsTCZ_Pre(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

#FGSEA_hl_Pre <- FGSEA_hl_Pre %>% mutate(Database = "Hallmark")
#FGSEA_wk_Pre <- FGSEA_wk_Pre %>% mutate(Database = "Wiki")
#FGSEA_kg_Pre <- FGSEA_kg_Pre %>% mutate(Database = "KEGG")

#combine_Pre <- bind_rows(FGSEA_hl_Pre, FGSEA_wk_Pre, FGSEA_kg_Pre)
#combine_Pre <- combine_Pre %>% select(Database, everything())

#wb <- loadWorkbook("out_fgsea_study1.xlsx")
#addWorksheet(wb, "TNFvsTCZ_Pre")
#writeData(wb, "TNFvsTCZ_Pre", combine_Pre)
#saveWorkbook(wb, "out_fgsea_Study1.xlsx", overwrite = T)

##--------------------------------------Barplots- FGSEA - TNFvsTCZ_Pre------------------------------------------
top_hl_Pre <- FGSEA_hl_Pre[pathway %in% topPathways_hl_Pre][order(pval)] #pull the top pathways from hallmark db

#Barplot for hallmark pathway - TNFvsTCZ_Pre

#jpeg("FGSEA_Hallmark_TNFvsTCZ_Pre.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_hl_Pre) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Hallmark Pathway - TNFvsTCZ_Pre")
#dev.off()


## Wiki pathways
top_wk_Pre <- FGSEA_wk_Pre[pathway %in% topPathways_wk_Pre][order(pval)] #pull the top pathways from wiki db

#Barplot for wiki pathway - TNFvsTCZ_Pre

#jpeg("FGSEA_Wiki_TNFvsTCZ_Pre.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_wk_Pre) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Wiki Pathway - TNFvsTCZ_Pre")
#dev.off()

## KEGG pathways
top_kg_Pre <- FGSEA_kg_Pre[pathway %in% topPathways_kg_Pre][order(pval)] #pull the top pathways from KEGG db

#Barplot for KEGG pathway - TNFvsTCZ_Pre

#jpeg("FGSEA_KEGG_TNFvsTCZ_Pre.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_kg_Pre) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="KEGG Pathway - TNFvsTCZ_Pre")
#dev.off()



##--------------------------------------Gene Set Enrichment Analysis - TNFvsTCZ_Post---------------------------
Hallmark <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt") # from MSigDB
Wiki <- gmtPathways("c2.cp.wikipathways.v2024.1.Hs.symbols.gmt") #from MSigDB
KEGG <- gmtPathways("c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt") #from MSigDB

# Prepare ranked gene list for Pre vs Post
ranked_genes_Post <- TNFvsTCZ_post_out$logFC
names(ranked_genes_Post) <- TNFvsTCZ_post_out$gene_name

# Sort genes in decreasing order of logFC
ranked_genes_Post <- sort(ranked_genes_Post, decreasing = TRUE)

# View the ranked genes
head(ranked_genes_Post)

# Run fgsea with Hallmark pathways 
FGSEA_hl_Post <- fgsea(pathways = Hallmark , # List of gene sets to check
                      stats = ranked_genes_Post,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_hl_Post$leadingEdge <- sapply(FGSEA_hl_Post$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_hl_Post, "Fgsea_Hallmark_Post.xlsx")

#select top pathways and create a table
topPathwaysUp_hl_Post <- FGSEA_hl_Post[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 hallmark pathways upregulated at pval 0.01
topPathwaysDown_hl_Post <- FGSEA_hl_Post[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 hallmark pathways downregulated at pval 0.01
topPathways_hl_Post <- c(topPathwaysUp_hl_Post, rev(topPathwaysDown_hl_Post))
plotGseaTable(Hallmark[topPathways_hl_Post], stats = ranked_genes_Post, fgseaRes = FGSEA_hl_Post, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_hl_Post <- collapsePathways(FGSEA_hl_Post[order(pval)][padj < 0.01], Hallmark, ranked_genes_Post, gseaParam = 1)

# plot the most significantly enriched pathway
plotEnrichment(Hallmark[[head(FGSEA_hl_Post[order(pval), ], 1)$pathway]],
               ranked_genes_Post) + 
  labs(title = head(FGSEA_hl_Post[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_hl_Post = Hallmark[topPathways_hl_Post]

# generate data for input into pathway images
df_hl_Post = gseaCurve(ranked_genes_Post, setlist_hl_Post, FGSEA_hl_Post)

#jpeg("HallmarkPathway_TNFvsTCZ_Post.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_hl_Post, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "HallmarkPathway_TNFvsTCZ_Post(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()


## Wiki pathways
# Run fgsea with Wiki pathways 
FGSEA_wk_Post <- fgsea(pathways = Wiki , # List of gene sets to check
                      stats = ranked_genes_Post,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_wk_Post$leadingEdge <- sapply(FGSEA_wk_Post$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_wk_Post, "Fgsea_Wiki_Post.xlsx")

#select top pathways and create a table
topPathwaysUp_wk_Post <- FGSEA_wk_Post[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 wiki pathways upregulated at pval 0.01
topPathwaysDown_wk_Post <- FGSEA_wk_Post[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 wiki pathways downregulated at pval 0.01
topPathways_wk_Post <- c(topPathwaysUp_wk_Post, rev(topPathwaysDown_wk_Post))
plotGseaTable(Wiki[topPathways_wk_Post], stats = ranked_genes_Post, fgseaRes = FGSEA_wk_Post, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_wk_Post <- collapsePathways(FGSEA_wk_Post[order(pval)][padj < 0.01], Wiki, ranked_genes_Post, gseaParam = 1)
mainPathways_wk_Post <- FGSEA_wk_Post[pathway %in% collapsedPathways_wk_Post$mainPathways][order(-NES), pathway]
plotGseaTable(Wiki[mainPathways_wk_Post], ranked_genes_Post, FGSEA_wk_Post, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(Wiki[[head(FGSEA_wk_Post[order(pval), ], 1)$pathway]],
               ranked_genes_Post) + 
  labs(title = head(FGSEA_wk_Post[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_wk_Post = Wiki[topPathways_wk_Post]

# generate data for input into pathway images
df_wk_Post = gseaCurve(ranked_genes_Post, setlist_wk_Post, FGSEA_wk_Post)

#jpeg("WikiPathway_TNFvsTCZ_Post.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_wk_Post, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "WikiPathway_TNFvsTCZ_Post(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

## KEGG pathway
# Run fgsea with KEGG pathways 
FGSEA_kg_Post <- fgsea(pathways = KEGG , # List of gene sets to check
                      stats = ranked_genes_Post,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_kg_Post$leadingEdge <- sapply(FGSEA_kg_Post$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_kg_Post, "Fgsea_KEGG_Post.xlsx")

#select top pathways and create a table
topPathwaysUp_kg_Post <- FGSEA_kg_Post[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 kegg pathways upregulated at pval 0.01
topPathwaysDown_kg_Post <- FGSEA_kg_Post[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 kegg pathways downregulated at pval 0.01
topPathways_kg_Post <- c(topPathwaysUp_kg_Post, rev(topPathwaysDown_kg_Post))
plotGseaTable(KEGG[topPathways_kg_Post], stats = ranked_genes_Post, fgseaRes = FGSEA_kg_Post, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_kg_Post <- collapsePathways(FGSEA_kg_Post[order(pval)][padj < 0.01], KEGG, ranked_genes_Post, gseaParam = 1)
mainPathways_kg_Post <- FGSEA_kg_Post[pathway %in% collapsedPathways_kg_Post$mainPathways][order(-NES), pathway]
plotGseaTable(KEGG[mainPathways_kg_Post], ranked_genes_Post, FGSEA_kg_Post, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(KEGG[[head(FGSEA_kg_Post[order(pval), ], 1)$pathway]],
               ranked_genes_Post) + 
  labs(title = head(FGSEA_kg_Post[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_kg_Post = KEGG[topPathways_kg_Post]

# generate data for input into pathway images
df_kg_Post = gseaCurve(ranked_genes_Post, setlist_kg_Post, FGSEA_kg_Post)

#jpeg("KEGGPathway_TNFvsTCZ_Post.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_kg_Post, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "KEGGPathway_TNFvsTCZ_Post(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

#FGSEA_hl_Post <- FGSEA_hl_Post %>% mutate(Database = "Hallmark")
#FGSEA_wk_Post <- FGSEA_wk_Post %>% mutate(Database = "Wiki")
#FGSEA_kg_Post <- FGSEA_kg_Post %>% mutate(Database = "KEGG")

#combine_Post <- bind_rows(FGSEA_hl_Post, FGSEA_wk_Post, FGSEA_kg_Post)
#combine_Post <- combine_Post %>% select(Database, everything())

#wb <- loadWorkbook("out_fgsea_study1.xlsx")
#addWorksheet(wb, "TNFvsTCZ_Post")
#writeData(wb, "TNFvsTCZ_Post", combine_Post)
#saveWorkbook(wb, "out_fgsea_Study1.xlsx", overwrite = T)

##--------------------------------------Barplots- FGSEA - TNFvsTCZ_Post------------------------------------------
top_hl_Post <- FGSEA_hl_Post[pathway %in% topPathways_hl_Post][order(pval)] #pull the top pathways from hallmark db

#Barplot for hallmark pathway - TNFvsTCZ_Post

#jpeg("FGSEA_Hallmark_TNFvsTCZ_Post.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_hl_Post) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Hallmark Pathway - TNFvsTCZ_Post")
#dev.off()


## Wiki pathways
top_wk_Post <- FGSEA_wk_Post[pathway %in% topPathways_wk_Post][order(pval)] #pull the top pathways from wiki db

#Barplot for wiki pathway - TNFvsTCZ_Post

#jpeg("FGSEA_Wiki_TNFvsTCZ_Post.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_wk_Post) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Wiki Pathway - TNFvsTCZ_Post")
#dev.off()

## KEGG pathways
top_kg_Post <- FGSEA_kg_Post[pathway %in% topPathways_kg_Post][order(pval)] #pull the top pathways from KEGG db

#Barplot for KEGG pathway - TNFvsTCZ_Post

#jpeg("FGSEA_KEGG_TNFvsTCZ_Post.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_kg_Post) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="KEGG Pathway - TNFvsTCZ_Post")
#dev.off()


##--------------------------------------Gene Set Enrichment Analysis - TNFvsTCZ_Intxn---------------------------
Hallmark <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt") # from MSigDB
Wiki <- gmtPathways("c2.cp.wikipathways.v2024.1.Hs.symbols.gmt") #from MSigDB
KEGG <- gmtPathways("c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt") #from MSigDB

# Prepare ranked gene list for Pre vs Post
ranked_genes_Intxn <- TNFvsTCZ_Intxn_out$logFC
names(ranked_genes_Intxn) <- TNFvsTCZ_Intxn_out$gene_name

# Sort genes in decreasing order of logFC
ranked_genes_Intxn <- sort(ranked_genes_Intxn, decreasing = TRUE)

# View the ranked genes
head(ranked_genes_Intxn)

# Run fgsea with Hallmark pathways 
FGSEA_hl_Intxn <- fgsea(pathways = Hallmark , # List of gene sets to check
                       stats = ranked_genes_Intxn,
                       eps = 0,
                       scoreType = 'std', 
                       minSize = 15,
                       maxSize = 250)

#FGSEA_hl_Intxn$leadingEdge <- sapply(FGSEA_hl_Intxn$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_hl_Intxn, "Fgsea_Hallmark_Intxn.xlsx")

#select top pathways and create a table
topPathwaysUp_hl_Intxn <- FGSEA_hl_Intxn[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 hallmark pathways upregulated at pval 0.01
topPathwaysDown_hl_Intxn <- FGSEA_hl_Intxn[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 hallmark pathways downregulated at pval 0.01
topPathways_hl_Intxn <- c(topPathwaysUp_hl_Intxn, rev(topPathwaysDown_hl_Intxn))
plotGseaTable(Hallmark[topPathways_hl_Intxn], stats = ranked_genes_Intxn, fgseaRes = FGSEA_hl_Intxn, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_hl_Intxn <- collapsePathways(FGSEA_hl_Intxn[order(pval)][padj < 0.01], Hallmark, ranked_genes_Intxn, gseaParam = 1)

# plot the most significantly enriched pathway
plotEnrichment(Hallmark[[head(FGSEA_hl_Intxn[order(pval), ], 1)$pathway]],
               ranked_genes_Intxn) + 
  labs(title = head(FGSEA_hl_Intxn[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_hl_Intxn = Hallmark[topPathways_hl_Intxn]

# generate data for input into pathway images
df_hl_Intxn = gseaCurve(ranked_genes_Intxn, setlist_hl_Intxn, FGSEA_hl_Intxn)

#jpeg("HallmarkPathway_TNFvsTCZ_Intxn.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_hl_Intxn, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "HallmarkPathway_TNFvsTCZ_Intxn(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()


## Wiki pathways
# Run fgsea with Wiki pathways 
FGSEA_wk_Intxn <- fgsea(pathways = Wiki , # List of gene sets to check
                       stats = ranked_genes_Intxn,
                       eps = 0,
                       scoreType = 'std', 
                       minSize = 15,
                       maxSize = 250)

#FGSEA_wk_Intxn$leadingEdge <- sapply(FGSEA_wk_Intxn$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_wk_Intxn, "Fgsea_Wiki_Intxn.xlsx")

#select top pathways and create a table
topPathwaysUp_wk_Intxn <- FGSEA_wk_Intxn[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 wiki pathways upregulated at pval 0.01
topPathwaysDown_wk_Intxn <- FGSEA_wk_Intxn[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 wiki pathways downregulated at pval 0.01
topPathways_wk_Intxn <- c(topPathwaysUp_wk_Intxn, rev(topPathwaysDown_wk_Intxn))
plotGseaTable(Wiki[topPathways_wk_Intxn], stats = ranked_genes_Intxn, fgseaRes = FGSEA_wk_Intxn, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_wk_Intxn <- collapsePathways(FGSEA_wk_Intxn[order(pval)][padj < 0.01], Wiki, ranked_genes_Intxn, gseaParam = 1)
mainPathways_wk_Intxn <- FGSEA_wk_Intxn[pathway %in% collapsedPathways_wk_Intxn$mainPathways][order(-NES), pathway]
plotGseaTable(Wiki[mainPathways_wk_Intxn], ranked_genes_Intxn, FGSEA_wk_Intxn, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(Wiki[[head(FGSEA_wk_Intxn[order(pval), ], 1)$pathway]],
               ranked_genes_Intxn) + 
  labs(title = head(FGSEA_wk_Intxn[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_wk_Intxn = Wiki[topPathways_wk_Intxn]

# generate data for input into pathway images
df_wk_Intxn = gseaCurve(ranked_genes_Intxn, setlist_wk_Intxn, FGSEA_wk_Intxn)

#jpeg("WikiPathway_TNFvsTCZ_Intxn.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_wk_Intxn, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "WikiPathway_TNFvsTCZ_Intxn(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

## KEGG pathway
# Run fgsea with KEGG pathways 
FGSEA_kg_Intxn <- fgsea(pathways = KEGG , # List of gene sets to check
                       stats = ranked_genes_Intxn,
                       eps = 0,
                       scoreType = 'std', 
                       minSize = 15,
                       maxSize = 250)

#FGSEA_kg_Intxn$leadingEdge <- sapply(FGSEA_kg_Intxn$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_kg_Intxn, "Fgsea_KEGG_Intxn.xlsx")

#select top pathways and create a table
topPathwaysUp_kg_Intxn <- FGSEA_kg_Intxn[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 kegg pathways upregulated at pval 0.01
topPathwaysDown_kg_Intxn <- FGSEA_kg_Intxn[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 kegg pathways downregulated at pval 0.01
topPathways_kg_Intxn <- c(topPathwaysUp_kg_Intxn, rev(topPathwaysDown_kg_Intxn))
plotGseaTable(KEGG[topPathways_kg_Intxn], stats = ranked_genes_Intxn, fgseaRes = FGSEA_kg_Intxn, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_kg_Intxn <- collapsePathways(FGSEA_kg_Intxn[order(pval)][padj < 0.01], KEGG, ranked_genes_Intxn, gseaParam = 1)
mainPathways_kg_Intxn <- FGSEA_kg_Intxn[pathway %in% collapsedPathways_kg_Intxn$mainPathways][order(-NES), pathway]
plotGseaTable(KEGG[mainPathways_kg_Intxn], ranked_genes_Intxn, FGSEA_kg_Intxn, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(KEGG[[head(FGSEA_kg_Intxn[order(pval), ], 1)$pathway]],
               ranked_genes_Intxn) + 
  labs(title = head(FGSEA_kg_Intxn[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_kg_Intxn = KEGG[topPathways_kg_Intxn]

# generate data for input into pathway images
df_kg_Intxn = gseaCurve(ranked_genes_Intxn, setlist_kg_Intxn, FGSEA_kg_Intxn)

#jpeg("KEGGPathway_TNFvsTCZ_Intxn.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_kg_Intxn, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "KEGGPathway_TNFvsTCZ_Intxn(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

#FGSEA_hl_Intxn <- FGSEA_hl_Intxn %>% mutate(Database = "Hallmark")
#FGSEA_wk_Intxn <- FGSEA_wk_Intxn %>% mutate(Database = "Wiki")
#FGSEA_kg_Intxn <- FGSEA_kg_Intxn %>% mutate(Database = "KEGG")

#combine_Intxn <- bind_rows(FGSEA_hl_Intxn, FGSEA_wk_Intxn, FGSEA_kg_Intxn)
#combine_Intxn <- combine_Intxn %>% select(Database, everything())

#wb <- loadWorkbook("out_fgsea_study1.xlsx")
#addWorksheet(wb, "TNFvsTCZ_Intxn")
#writeData(wb, "TNFvsTCZ_Intxn", combine_Intxn)
#saveWorkbook(wb, "out_fgsea_Study1.xlsx", overwrite = T)

##--------------------------------------Barplots- FGSEA - TNFvsTCZ_Intxn------------------------------------------
top_hl_Intxn <- FGSEA_hl_Intxn[pathway %in% topPathways_hl_Intxn][order(pval)] #pull the top pathways from hallmark db

#Barplot for hallmark pathway - TNFvsTCZ_Intxn

#jpeg("FGSEA_Hallmark_TNFvsTCZ_Intxn.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_hl_Intxn) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Hallmark Pathway - TNFvsTCZ_Intxn")
#dev.off()


## Wiki pathways
top_wk_Intxn <- FGSEA_wk_Intxn[pathway %in% topPathways_wk_Intxn][order(pval)] #pull the top pathways from wiki db

#Barplot for wiki pathway - TNFvsTCZ_Intxn

#jpeg("FGSEA_Wiki_TNFvsTCZ_Intxn.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_wk_Intxn) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Wiki Pathway - TNFvsTCZ_Intxn")
#dev.off()

## KEGG pathways
top_kg_Intxn <- FGSEA_kg_Intxn[pathway %in% topPathways_kg_Intxn][order(pval)] #pull the top pathways from KEGG db

#Barplot for KEGG pathway - TNFvsTCZ_Intxn

#jpeg("FGSEA_KEGG_TNFvsTCZ_Intxn.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_kg_Intxn) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="KEGG Pathway - TNFvsTCZ_Intxn")
#dev.off()



##-----------------------------------------------------paired analysis - TNF group---------------
## TNF data 
TNF_raw_counts <- study1_raw_counts_out[,1:7]
#colnames(TNF_raw_counts)[1:7] <- samples.time[1:8]

TNF_DGE <- DGEList(counts = TNF_raw_counts) # convert the raw counts into DGE list

## TMM normalization
TNF_DGE <- calcNormFactors(TNF_DGE)

## log2 transformation of count data and prior count = 3
TNF_logCPM <- cpm(TNF_DGE$counts, log = T, prior.count = 3) # logcpm values without gene names

TNF_logCPM_file <- data.frame(Gene = TNF_DGE$genes$gene_name, TNF_logCPM) #logcpm values with gene names
#write.xlsx(TNF_logCPM_file, "TNF_logCPM.xlsx")

TNF_DGE$samples$group <- study1_groups_out$Groups[1:6] # add the TNF group to the DGE list

# Filter the lowly expressed genes
TNF_keep <- rowSums(cpm(TNF_DGE) >1) >= 3
TNF_DGE <- TNF_DGE[TNF_keep,,keep.lib.sizes = F]
length(rownames(TNF_DGE$counts)) #12629 genes remained after filtering

##Differential gene expression analysis - TNF data
TNF_subject <- factor(study1_groups$Patient_Combined[c(1,2,5:8)])
TNF_subject <- factor(TNF_subject, levels = unique(TNF_subject)) #get the subject data

## Time 
TNF_Time <- factor(study1_groups$Time[c(1,2,5:8)])
TNF_Time <- factor(TNF_Time, levels = unique(TNF_Time))

TNF_paired_data <- data.frame(Samplenames = colnames(TNF_DGE$counts), Time = TNF_Time, Subject = TNF_subject)

## create a design matrix for the time comparison
TNF_design <- model.matrix(~0+Time+Subject, data = TNF_paired_data)
colnames(TNF_design) <- gsub("Time", "", colnames(TNF_design))
colnames(TNF_design) <- gsub("Subject", "", colnames(TNF_design))
rownames(TNF_design) <- TNF_paired_data$Samplenames

## contrast matrix for time comparison
TNF_contr.mtrx <- makeContrasts(
  Time = post - pre,
  levels = TNF_design)
TNF_contr.mtrx

## run voom
TNF_v <- voom(TNF_DGE, TNF_design, plot = T)
TNF_fit <- lmFit(TNF_v, TNF_design)
TNF_fit <- contrasts.fit(TNF_fit, TNF_contr.mtrx)
TNF_efit <- eBayes(TNF_fit)

TNF_PostvsPre <- topTable(TNF_efit, n = Inf, adjust.method = "fdr", sort.by = "P") # no fdr threshold was given to this data - contains all genes
TNF_PostvsPre_FDR_0.05 <- topTable(TNF_efit, n = Inf, adjust.method = "fdr", sort.by = "P", p.value = 0.05) #found 2286 genes at FDR 0.05
TNF_PostvsPre_FDR_0.01 <- topTable(TNF_efit, n = Inf, adjust.method = "fdr", sort.by = "P", p.value = 0.01) #found 1120 genes at FDR 0.01

#write.xlsx(TNF_PostvsPre, "out_limma_tnf.xlsx")

## Volcano plot
#jpeg("VolcanoPlot_TNF_PostvsPre.jpeg", width = 2500, height = 1800, res = 300)
ggplot(TNF_PostvsPre, aes(x= logFC ,
                              y=-log10(adj.P.Val), 
                              color= -log10(adj.P.Val)))+
  geom_point(size=2)+
  geom_point(data=subset(TNF_PostvsPre %>% slice_min(adj.P.Val, n=20)),
             aes(x= logFC, 
                 y=-log10(adj.P.Val), 
                 color= -log10(adj.P.Val)))+
  scale_color_viridis(option="rocket", direction=1)+
  geom_text_repel(data=subset(TNF_PostvsPre %>% slice_min(adj.P.Val, n=20)),
                  aes(x=logFC, 
                      y=-log10(adj.P.Val), 
                      color= -log10(adj.P.Val),label = gene_name), 
                  color = "black", size=4, box.padding = unit(0.2, "lines"),segment.color="gray25",
                  segment.size=0.25, point.padding = unit(0.2, "lines"), max.overlaps = Inf)+
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(1,-1),linetype="dashed")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text (size = 10, face ="bold", colour = "black"))+
  theme(axis.title = element_text (size = 18, face = "bold", colour = "black")) +
  theme(strip.text = element_text(size=10, face="bold", color = "black"))+
  ggtitle("volcano_Plot_Top20genes_TNF_PostvsPre")+
  theme(plot.title=element_text(face="bold"))+
  xlab("log2FC (Post vs. Pre)") +
  ylab("-log10 FDR (Post vs. Pre)")+
  theme(legend.position="bottom") #to maximize graph space by removing legends
#dev.off()

## Heatmaps for paired analysis - TNF data
TNF_PostvsPre_filtered <- TNF_PostvsPre[TNF_PostvsPre$adj.P.Val < 0.05 & abs(TNF_PostvsPre$logFC) > 1, ]
TNF.i.1 <- TNF_efit$genes$gene_name[which(TNF_efit$genes$gene_name %in% TNF_PostvsPre_filtered$gene_name)] 

data_TNF_PostvsPre <- TNF_logCPM_file[TNF_logCPM_file$Gene %in% TNF.i.1, -1]
rownames(data_TNF_PostvsPre) <- TNF.i.1

TNF.samples.time <- paste(colnames(data_TNF_PostvsPre), TNF_Time, sep = "_")
colnames(data_TNF_PostvsPre) <- TNF.samples.time

TNF_pre_samples <- grep("pre", colnames(data_TNF_PostvsPre), value = TRUE)  # Extract pre samples
TNF_post_samples <- grep("post", colnames(data_TNF_PostvsPre), value = T) #extract post samples

# Reorder the columns such that pre and post samples are grouped
TNF_ordered_samples <- c(TNF_pre_samples, TNF_post_samples)

#color palette
palette <- colorRampPalette(c("blue", "white", "red"))(100)

data_TNF_PostvsPre <- data_TNF_PostvsPre[, TNF_ordered_samples]

# create an annotation dataframe for the samples
annot_col_TNF_PrevsPost <- data.frame(Time = substr(colnames(data_TNF_PostvsPre),6,9))
rownames(annot_col_TNF_PrevsPost) <- colnames(data_TNF_PostvsPre)


annot_colors_TNF_PrevsPost <- list(
  Time = c(pre = "#1B9E77", post = "#D95F02")
)

data_TNF_PostvsPre <- data_TNF_PostvsPre %>% slice(1:100)

#jpeg("Heatmap_TNF_PostvsPre.jpeg", width = 3000, height = 3500, res = 300)
#generate heatmap
pheatmap(data_TNF_PostvsPre, show_rownames = T, scale = "row", cluster_cols= F, cluster_rows = T, color = palette, border_color = NA,
         fontsize_row = 8 , fontsize_col = 10,annotation_col = annot_col_TNF_PrevsPost, annotation_colors = annot_colors_TNF_PrevsPost, 
         annotation_legend = T, main = "Topgenes_TNF_PrevsPost_FDR_0.05")
#dev.off()

##---------------------------------------------------------Overrepresentation analysis - TNF group----------------------------------------
#TNF_genes_FDR_0.05 <- data.frame(gene_name = TNF_PostvsPre_FDR_0.05$gene_name)
#write.xlsx(TNF_genes_FDR_0.05, "TNF_genes.xlsx")

# Load the pathway tables extracted from enrichr
Hallmark_TNF <- read.delim("Hallmark_TNF.txt", sep = '\t')
Wiki_TNF <- read.delim("WikiPathways_TNF.txt", sep = '\t')
KEGG_TNF <- read.delim("KEGG_TNF.txt", sep = '\t')

#write.xlsx(Hallmark_TNF, "out_enrichr_TNF.xlsx")
#wb_tnf <- loadWorkbook("out_enrichr_TNF.xlsx")
#addWorksheet(wb_tnf, "KEGGPathways")
#writeData(wb_tnf, "KEGGPathways", KEGG_TNF)
#saveWorkbook(wb_tnf, "out_enrichr_TNF.xlsx", overwrite = T)

# filter the hallmark, wiki and kegg pathways based on adj.p.val less than 0.05
Hallmark_TNF <- Hallmark_TNF %>% filter(Adjusted.P.value < 0.05)
Wiki_TNF <- Wiki_TNF %>% filter(Adjusted.P.value < 0.05)
KEGG_TNF <- KEGG_TNF %>% filter(Adjusted.P.value < 0.05)

#Plot the filtered hallmark pathways
#jpeg("HallmarkPathways_TNF_PostvsPre.jpeg", width = 2500, height = 1800, res = 300)
ggplot(Hallmark_TNF) +
  geom_col(aes(
    x = reorder(Term, Odds.Ratio),
    y = Odds.Ratio,
    fill = -log10(Adjusted.P.value)
  ), color = "black") +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  coord_flip() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  labs(y = "Odds Ratio", 
       x = "Pathways", 
       title = "HallmarkPathway - TNF_PostvsPre (Adj.P.Val < 0.05)",
       fill = "-log10(Adjusted.P.value)")  
#dev.off()

#Plot the filtered wiki pathways
#jpeg("WikiPathways_TNF_PostvsPre.jpeg", width = 2500, height = 1800, res = 300)
ggplot(Wiki_TNF %>% slice_min(Adjusted.P.value, n = 50)) +
  geom_col(aes(
    x = reorder(Term, Odds.Ratio),
    y = Odds.Ratio,
    fill = -log10(Adjusted.P.value)
  ), color = "black") +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  coord_flip() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  labs(y = "Odds Ratio", 
       x = "Pathways", 
       title = "WikiPathway - TNF_PostvsPre (Adj.P.Val < 0.05)",
       fill = "-log10(Adjusted.P.value)")
#dev.off()

#Plot the filtered kegg pathways
#jpeg("KEGGPathways_TNF_PostvsPre.jpeg", width = 3500, height = 2800, res = 300)
ggplot(KEGG_TNF) +
  geom_col(aes(
    x = reorder(Term, Odds.Ratio),
    y = Odds.Ratio,
    fill = -log10(Adjusted.P.value)
  ), color = "black") +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  coord_flip() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  labs(y = "Odds Ratio", 
       x = "Pathways", 
       title = "KEGGPathway - TNF_PostvsPre (Adj.P.Val < 0.05)",
       fill = "-log10(Adjusted.P.value)")  
#dev.off()


##------------------------------------------------Gene Set Enrichement Analysis - TNF---------------------------------
Hallmark <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt") # from MSigDB
Wiki <- gmtPathways("c2.cp.wikipathways.v2024.1.Hs.symbols.gmt") #from MSigDB
KEGG <- gmtPathways("c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt") #from MSigDB

# prepare ranked gene list for Pre vs Post
ranked_genes_TNF <- TNF_PostvsPre$logFC
names(ranked_genes_TNF) <- TNF_PostvsPre$gene_name

# Sort genes in decreasing order of logFC
ranked_genes_TNF <- sort(ranked_genes_TNF, decreasing = TRUE)

# View the ranked genes
head(ranked_genes_TNF)

# Run fgsea with Hallmark pathways 
FGSEA_hl_TNF <- fgsea(pathways = Hallmark , # List of gene sets to check
                       stats = ranked_genes_TNF,
                       eps = 0,
                       scoreType = 'std', 
                       minSize = 15,
                       maxSize = 250)

#FGSEA_hl_TNF$leadingEdge <- sapply(FGSEA_hl_TNF$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_hl_TNF, "Fgsea_Hallmark_TNF.xlsx")

#select top pathways and create a table
topPathwaysUp_hl_TNF <- FGSEA_hl_TNF[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 hallmark pathways upregulated at pval 0.01
topPathwaysDown_hl_TNF <- FGSEA_hl_TNF[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 hallmark pathways downregulated at pval 0.01
topPathways_hl_TNF <- c(topPathwaysUp_hl_TNF, rev(topPathwaysDown_hl_TNF))
plotGseaTable(Hallmark[topPathways_hl_TNF], stats = ranked_genes_TNF, fgseaRes = FGSEA_hl_TNF, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_hl_TNF <- collapsePathways(FGSEA_hl_TNF[order(pval)][padj < 0.01], Hallmark, ranked_genes_TNF, gseaParam = 1)

# plot the most significantly enriched pathway
plotEnrichment(Hallmark[[head(FGSEA_hl_TNF[order(pval), ], 1)$pathway]],
               ranked_genes_TNF) + 
  labs(title = head(FGSEA_hl_TNF[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_hl_TNF = Hallmark[topPathways_hl_TNF]

# generate data for input into pathway images
df_hl_TNF = gseaCurve(ranked_genes_TNF, setlist_hl_TNF, FGSEA_hl_TNF)

#jpeg("HallmarkPathway_TNF.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_hl_TNF, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "HallmarkPathway_TNF_PostvsPre(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()


## Wiki pathways
# Run fgsea with Wiki pathways 
FGSEA_wk_TNF <- fgsea(pathways = Wiki , # List of gene sets to check
                      stats = ranked_genes_TNF,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_wk_TNF$leadingEdge <- sapply(FGSEA_wk_TNF$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_wk_TNF, "Fgsea_Wiki_TNF.xlsx")

#select top pathways and create a table
topPathwaysUp_wk_TNF <- FGSEA_wk_TNF[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 wiki pathways upregulated at pval 0.01
topPathwaysDown_wk_TNF <- FGSEA_wk_TNF[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 wiki pathways downregulated at pval 0.01
topPathways_wk_TNF <- c(topPathwaysUp_wk_TNF, rev(topPathwaysDown_wk_TNF))
plotGseaTable(Wiki[topPathways_wk_TNF], stats = ranked_genes_TNF, fgseaRes = FGSEA_wk_TNF, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_wk_TNF <- collapsePathways(FGSEA_wk_TNF[order(pval)][padj < 0.01], Wiki, ranked_genes_TNF, gseaParam = 1)
mainPathways_wk_TNF <- FGSEA_wk_TNF[pathway %in% collapsedPathways_wk_TNF$mainPathways][order(-NES), pathway]
plotGseaTable(Wiki[mainPathways_wk_TNF], ranked_genes_TNF, FGSEA_wk_TNF, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(Wiki[[head(FGSEA_wk_TNF[order(pval), ], 1)$pathway]],
               ranked_genes_TNF) + 
  labs(title = head(FGSEA_wk_TNF[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_wk_TNF = Wiki[topPathways_wk_TNF]

# generate data for input into pathway images
df_wk_TNF = gseaCurve(ranked_genes_TNF, setlist_wk_TNF, FGSEA_wk_TNF)

#jpeg("WikiPathway_TNF.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_wk_TNF, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "WikiPathway_TNF_PostvsPre(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

## KEGG pathway
# Run fgsea with KEGG pathways 
FGSEA_kg_TNF <- fgsea(pathways = KEGG , # List of gene sets to check
                      stats = ranked_genes_TNF,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_kg_TNF$leadingEdge <- sapply(FGSEA_kg_TNF$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_kg_TNF, "Fgsea_KEGG_TNF.xlsx")

#select top pathways and create a table
topPathwaysUp_kg_TNF <- FGSEA_kg_TNF[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 kegg pathways upregulated at pval 0.01
topPathwaysDown_kg_TNF <- FGSEA_kg_TNF[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 kegg pathways downregulated at pval 0.01
topPathways_kg_TNF <- c(topPathwaysUp_kg_TNF, rev(topPathwaysDown_kg_TNF))
plotGseaTable(KEGG[topPathways_kg_TNF], stats = ranked_genes_TNF, fgseaRes = FGSEA_kg_TNF, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_kg_TNF <- collapsePathways(FGSEA_kg_TNF[order(pval)][padj < 0.01], KEGG, ranked_genes_TNF, gseaParam = 1)
mainPathways_kg_TNF <- FGSEA_kg_TNF[pathway %in% collapsedPathways_kg_TNF$mainPathways][order(-NES), pathway]
plotGseaTable(KEGG[mainPathways_kg_TNF], ranked_genes_TNF, FGSEA_kg_TNF, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(KEGG[[head(FGSEA_kg_TNF[order(pval), ], 1)$pathway]],
               ranked_genes_TNF) + 
  labs(title = head(FGSEA_kg_TNF[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_kg_TNF = KEGG[topPathways_kg_TNF]

# generate data for input into pathway images
df_kg_TNF = gseaCurve(ranked_genes_TNF, setlist_kg_TNF, FGSEA_kg_TNF)

#jpeg("KEGGPathway_TNF.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_kg_TNF, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "KEGGPathway_TNF_PostvsPre(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

##--------------------------------------Barplots- FGSEA - TNF------------------------------------------
top_hl_TNF <- FGSEA_hl_TNF[pathway %in% topPathways_hl_TNF][order(pval)] #pull the top pathways from hallmark db

#Barplot for hallmark pathway - TNF

#jpeg("FGSEA_Hallmark_TNF.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_hl_TNF) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Hallmark Pathway - TNF")
#dev.off()


## Wiki pathways
top_wk_TNF <- FGSEA_wk_TNF[pathway %in% topPathways_wk_TNF][order(pval)] #pull the top pathways from wiki db

#Barplot for wiki pathway - TNF

#jpeg("FGSEA_Wiki_TNF.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_wk_TNF) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Wiki Pathway - TNF")
#dev.off()

## KEGG pathways
top_kg_TNF <- FGSEA_kg_TNF[pathway %in% topPathways_kg_TNF][order(pval)] #pull the top pathways from KEGG db

#Barplot for KEGG pathway - TNF

#jpeg("FGSEA_KEGG_TNF.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_kg_TNF) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="KEGG Pathway - TNF")
#dev.off()


##--------------------------------------------------paired analysis - TCZ group------------------------

## TCZ data 
TCZ_raw_counts <- study1_raw_counts_out[,c(1,8:19)]
#colnames(TCZ_raw_counts)[2:15] <- samples.time[9:22]

TCZ_DGE <- DGEList(counts = TCZ_raw_counts) # convert the raw counts into DGE list

## TMM normalization
TCZ_DGE <- calcNormFactors(TCZ_DGE)

## log2 transformation of count data and prior count = 3
TCZ_logCPM <- cpm(TCZ_DGE$counts, log = T, prior.count = 3) # logcpm values without gene names

TCZ_logCPM_file <- data.frame(Gene = TCZ_DGE$genes$gene_name, TCZ_logCPM) #logcpm values with gene names
#write.xlsx(TCZ_logCPM_file, "TCZ_logCPM.xlsx")

TCZ_DGE$samples$group <- study1_groups_out$Groups[7:18] # add the TCZ group to the DGE list

# Filter the lowly expressed genes
TCZ_keep <- rowSums(cpm(TCZ_DGE) >1) >= 6
TCZ_DGE <- TCZ_DGE[TCZ_keep,,keep.lib.sizes = F]
length(rownames(TCZ_DGE$counts)) #12466 genes remained after filtering

##Differential gene expression analysis - TCZ data
TCZ_subject <- factor(study1_groups$Patient_Combined[9:20])
TCZ_subject <- factor(TCZ_subject, levels = unique(TCZ_subject)) #get the subject data

## Time 
TCZ_Time <- factor(study1_groups$Time[9:20])
TCZ_Time <- factor(TCZ_Time, levels = unique(TCZ_Time))

TCZ_paired_data <- data.frame(Samplenames = colnames(TCZ_DGE$counts), Time = TCZ_Time, Subject = TCZ_subject)

## create a design matrix for the time comparison
TCZ_design <- model.matrix(~0+Time+Subject, data = TCZ_paired_data)
colnames(TCZ_design) <- gsub("Time", "", colnames(TCZ_design))
colnames(TCZ_design) <- gsub("Subject", "", colnames(TCZ_design))
rownames(TCZ_design) <- TCZ_paired_data$Samplenames

## contrast matrix for time comparison
TCZ_contr.mtrx <- makeContrasts(
  Time = post - pre,
  levels = TCZ_design)
TCZ_contr.mtrx

## run voom
TCZ_v <- voom(TCZ_DGE, TCZ_design, plot = T)
TCZ_fit <- lmFit(TCZ_v, TCZ_design)
TCZ_fit <- contrasts.fit(TCZ_fit, TCZ_contr.mtrx)
TCZ_efit <- eBayes(TCZ_fit)

TCZ_PostvsPre <- topTable(TCZ_efit, n = Inf, adjust.method = "fdr", sort.by = "P") # no fdr threshold was given to this data - contains all genes
TCZ_PostvsPre_FDR_0.05 <- topTable(TCZ_efit, n = Inf, adjust.method = "fdr", sort.by = "P", p.value = 0.05) #found 2286 genes at FDR 0.05
TCZ_PostvsPre_FDR_0.01 <- topTable(TCZ_efit, n = Inf, adjust.method = "fdr", sort.by = "P", p.value = 0.01) #found 1120 genes at FDR 0.01

#write.xlsx(TCZ_PostvsPre, "out_limma_tcz.xlsx")

## Volcano plot
#jpeg("VolcanoPlot_TCZ_PostvsPre.jpeg", width = 2500, height = 1800, res = 300)
ggplot(TCZ_PostvsPre, aes(x= logFC ,
                          y=-log10(adj.P.Val), 
                          color= -log10(adj.P.Val)))+
  geom_point(size=2)+
  geom_point(data=subset(TCZ_PostvsPre %>% slice_min(adj.P.Val, n=20)),
             aes(x= logFC, 
                 y=-log10(adj.P.Val), 
                 color= -log10(adj.P.Val)))+
  scale_color_viridis(option="rocket", direction=1)+
  geom_text_repel(data=subset(TCZ_PostvsPre %>% slice_min(adj.P.Val, n=20)),
                  aes(x=logFC, 
                      y=-log10(adj.P.Val), 
                      color= -log10(adj.P.Val),label = gene_name), 
                  color = "black", size=4, box.padding = unit(0.2, "lines"),segment.color="gray25",
                  segment.size=0.25, point.padding = unit(0.2, "lines"), max.overlaps = Inf)+
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(1,-1),linetype="dashed")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text (size = 10, face ="bold", colour = "black"))+
  theme(axis.title = element_text (size = 18, face = "bold", colour = "black")) +
  theme(strip.text = element_text(size=10, face="bold", color = "black"))+
  ggtitle("volcano_Plot_Top20genes_TCZ_PostvsPre")+
  theme(plot.title=element_text(face="bold"))+
  xlab("log2FC (Post vs. Pre)") +
  ylab("-log10 FDR (Post vs. Pre)")+
  theme(legend.position="bottom") #to maximize graph space by removing legends
#dev.off()

## Heatmaps for paired analysis - TCZ data
TCZ_PostvsPre_filtered <- TCZ_PostvsPre[TCZ_PostvsPre$adj.P.Val < 0.05 & abs(TCZ_PostvsPre$logFC) > 1, ]
TCZ.i.1 <- TCZ_efit$genes$gene_name[which(TCZ_efit$genes$gene_name %in% TCZ_PostvsPre_filtered$gene_name)] 

data_TCZ_PostvsPre <- TCZ_logCPM_file[TCZ_logCPM_file$Gene %in% TCZ.i.1, -1]
rownames(data_TCZ_PostvsPre) <- TCZ.i.1

TNF.samples.time <- paste(colnames(data_TCZ_PostvsPre), TCZ_Time, sep = "_")
colnames(data_TCZ_PostvsPre) <- TNF.samples.time

TCZ_pre_samples <- grep("pre", colnames(data_TCZ_PostvsPre), value = TRUE)  # Extract pre samples
TCZ_post_samples <- grep("post", colnames(data_TCZ_PostvsPre), value = T) #extract post samples

# Reorder the columns such that pre and post samples are grouped
TCZ_ordered_samples <- c(TCZ_pre_samples, TCZ_post_samples)

#color palette
palette <- colorRampPalette(c("blue", "white", "red"))(100)

data_TCZ_PostvsPre <- data_TCZ_PostvsPre[, TCZ_ordered_samples]

# create an annotation dataframe for the samples
annot_col_TCZ_PrevsPost <- data.frame(Time = substr(colnames(data_TCZ_PostvsPre),6,9))
rownames(annot_col_TCZ_PrevsPost) <- colnames(data_TCZ_PostvsPre)


annot_colors_TCZ_PrevsPost <- list(
  Time = c(pre = "#1B9E77", post = "#D95F02")
)

data_TCZ_PostvsPre <- data_TCZ_PostvsPre %>% slice(1:100)

#jpeg("Heatmap_TCZ_PostvsPre.jpeg", width = 3000, height = 3500, res = 300)
#generate heatmap
pheatmap(data_TCZ_PostvsPre, show_rownames = T, scale = "row", cluster_cols= F, cluster_rows = T, color = palette, border_color = NA,
         fontsize_row = 8 , fontsize_col = 10,annotation_col = annot_col_TCZ_PrevsPost, annotation_colors = annot_colors_TCZ_PrevsPost, 
         annotation_legend = T, main = "Topgenes_TCZ_PrevsPost_FDR_0.05")
#dev.off()


##------------------------------------------------Gene Set Enrichement Analysis - TCZ---------------------------------
Hallmark <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt") # from MSigDB
Wiki <- gmtPathways("c2.cp.wikipathways.v2024.1.Hs.symbols.gmt") #from MSigDB
KEGG <- gmtPathways("c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt") #from MSigDB

# prepare ranked gene list for Pre vs Post
ranked_genes_TCZ <- TCZ_PostvsPre$logFC
names(ranked_genes_TCZ) <- TCZ_PostvsPre$gene_name

# Sort genes in decreasing order of logFC
ranked_genes_TCZ <- sort(ranked_genes_TCZ, decreasing = TRUE)

# View the ranked genes
head(ranked_genes_TCZ)

# Run fgsea with Hallmark pathways 
FGSEA_hl_TCZ <- fgsea(pathways = Hallmark , # List of gene sets to check
                      stats = ranked_genes_TCZ,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_hl_TCZ$leadingEdge <- sapply(FGSEA_hl_TCZ$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_hl_TCZ, "Fgsea_Hallmark_TCZ.xlsx")

#select top pathways and create a table
topPathwaysUp_hl_TCZ <- FGSEA_hl_TCZ[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 hallmark pathways upregulated at pval 0.01
topPathwaysDown_hl_TCZ <- FGSEA_hl_TCZ[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 hallmark pathways downregulated at pval 0.01
topPathways_hl_TCZ <- c(topPathwaysUp_hl_TCZ, rev(topPathwaysDown_hl_TCZ))
plotGseaTable(Hallmark[topPathways_hl_TCZ], stats = ranked_genes_TCZ, fgseaRes = FGSEA_hl_TCZ, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_hl_TCZ <- collapsePathways(FGSEA_hl_TCZ[order(pval)][padj < 0.01], Hallmark, ranked_genes_TCZ, gseaParam = 1)

# plot the most significantly enriched pathway
plotEnrichment(Hallmark[[head(FGSEA_hl_TCZ[order(pval), ], 1)$pathway]],
               ranked_genes_TCZ) + 
  labs(title = head(FGSEA_hl_TCZ[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_hl_TCZ = Hallmark[topPathways_hl_TCZ]

# generate data for input into pathway images
df_hl_TCZ = gseaCurve(ranked_genes_TCZ, setlist_hl_TCZ, FGSEA_hl_TCZ)

#jpeg("HallmarkPathway_TCZ.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_hl_TCZ, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "HallmarkPathway_TCZ_PostvsPre(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()


## Wiki pathways
# Run fgsea with Wiki pathways 
FGSEA_wk_TCZ <- fgsea(pathways = Wiki , # List of gene sets to check
                      stats = ranked_genes_TCZ,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_wk_TCZ$leadingEdge <- sapply(FGSEA_wk_TCZ$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_wk_TCZ, "Fgsea_Wiki_TCZ.xlsx")

#select top pathways and create a table
topPathwaysUp_wk_TCZ <- FGSEA_wk_TCZ[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 wiki pathways upregulated at pval 0.01
topPathwaysDown_wk_TCZ <- FGSEA_wk_TCZ[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 wiki pathways downregulated at pval 0.01
topPathways_wk_TCZ <- c(topPathwaysUp_wk_TCZ, rev(topPathwaysDown_wk_TCZ))
plotGseaTable(Wiki[topPathways_wk_TCZ], stats = ranked_genes_TCZ, fgseaRes = FGSEA_wk_TCZ, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_wk_TCZ <- collapsePathways(FGSEA_wk_TCZ[order(pval)][padj < 0.01], Wiki, ranked_genes_TCZ, gseaParam = 1)
mainPathways_wk_TCZ <- FGSEA_wk_TCZ[pathway %in% collapsedPathways_wk_TCZ$mainPathways][order(-NES), pathway]
plotGseaTable(Wiki[mainPathways_wk_TCZ], ranked_genes_TCZ, FGSEA_wk_TCZ, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(Wiki[[head(FGSEA_wk_TCZ[order(pval), ], 1)$pathway]],
               ranked_genes_TCZ) + 
  labs(title = head(FGSEA_wk_TCZ[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_wk_TCZ = Wiki[topPathways_wk_TCZ]

# generate data for input into pathway images
df_wk_TCZ = gseaCurve(ranked_genes_TCZ, setlist_wk_TCZ, FGSEA_wk_TCZ)

#jpeg("WikiPathway_TCZ.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_wk_TCZ, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "WikiPathway_TCZ_PostvsPre(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

## KEGG pathway
# Run fgsea with KEGG pathways 
FGSEA_kg_TCZ <- fgsea(pathways = KEGG , # List of gene sets to check
                      stats = ranked_genes_TCZ,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_kg_TCZ$leadingEdge <- sapply(FGSEA_kg_TCZ$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_kg_TCZ, "Fgsea_KEGG_TCZ.xlsx")

#select top pathways and create a table
topPathwaysUp_kg_TCZ <- FGSEA_kg_TCZ[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 kegg pathways upregulated at pval 0.01
topPathwaysDown_kg_TCZ <- FGSEA_kg_TCZ[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 kegg pathways downregulated at pval 0.01
topPathways_kg_TCZ <- c(topPathwaysUp_kg_TCZ, rev(topPathwaysDown_kg_TCZ))
plotGseaTable(KEGG[topPathways_kg_TCZ], stats = ranked_genes_TCZ, fgseaRes = FGSEA_kg_TCZ, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_kg_TCZ <- collapsePathways(FGSEA_kg_TCZ[order(pval)][padj < 0.01], KEGG, ranked_genes_TCZ, gseaParam = 1)
mainPathways_kg_TCZ <- FGSEA_kg_TCZ[pathway %in% collapsedPathways_kg_TCZ$mainPathways][order(-NES), pathway]
plotGseaTable(KEGG[mainPathways_kg_TCZ], ranked_genes_TCZ, FGSEA_kg_TCZ, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(KEGG[[head(FGSEA_kg_TCZ[order(pval), ], 1)$pathway]],
               ranked_genes_TCZ) + 
  labs(title = head(FGSEA_kg_TCZ[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_kg_TCZ = KEGG[topPathways_kg_TCZ]

# generate data for input into pathway images
df_kg_TCZ = gseaCurve(ranked_genes_TCZ, setlist_kg_TCZ, FGSEA_kg_TCZ)

#jpeg("KEGGPathway_TCZ.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_kg_TCZ, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "KEGGPathway_TCZ_PostvsPre(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

##--------------------------------------Barplots- FGSEA - TCZ------------------------------------------
top_hl_TCZ <- FGSEA_hl_TCZ[pathway %in% topPathways_hl_TCZ][order(pval)] #pull the top pathways from hallmark db

#Barplot for hallmark pathway - TCZ

#jpeg("FGSEA_Hallmark_TCZ.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_hl_TCZ) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Hallmark Pathway - TCZ")
#dev.off()


## Wiki pathways
top_wk_TCZ <- FGSEA_wk_TCZ[pathway %in% topPathways_wk_TCZ][order(pval)] #pull the top pathways from wiki db

#Barplot for wiki pathway - TCZ

#jpeg("FGSEA_Wiki_TCZ.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_wk_TCZ) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Wiki Pathway - TCZ")
#dev.off()

## KEGG pathways
top_kg_TCZ <- FGSEA_kg_TCZ[pathway %in% topPathways_kg_TCZ][order(pval)] #pull the top pathways from KEGG db

#Barplot for KEGG pathway - TCZ

#jpeg("FGSEA_KEGG_TCZ.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_kg_TCZ) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="KEGG Pathway - TCZ")
#dev.off()

##----------------------------------------------Overrepresentation analysis - TCZ --------------------------------------------------
#TCZ_genes_FDR_0.05 <- data.frame(gene_name = TCZ_PostvsPre_FDR_0.05$gene_name)
#write.xlsx(TCZ_genes_FDR_0.05, "TCZ_genes.xlsx")

Hallmark_TCZ <- read.delim("Hallmark_TCZ.txt", sep = '\t')
Wiki_TCZ <- read.delim("WikiPathways_TCZ.txt", sep = '\t')
KEGG_TCZ <- read.delim("KEGG_TCZ.txt", sep = '\t')

#write.xlsx(Hallmark_TCZ, "out_enrichr_TCZ.xlsx")
#wb_tcz <- loadWorkbook("out_enrichr_TCZ.xlsx")
#addWorksheet(wb_tcz, "KEGGPathways")
#writeData(wb_tcz, "KEGGPathways", KEGG_TCZ)
#saveWorkbook(wb_tcz, "out_enrichr_TCZ.xlsx", overwrite = T)


# filter the hallmark, wiki and kegg pathways based on adj.p.val less than 0.05
Hallmark_TCZ <- Hallmark_TCZ %>% filter(Adjusted.P.value < 0.05)
Wiki_TCZ <- Wiki_TCZ %>% filter(Adjusted.P.value < 0.05)
KEGG_TCZ <- KEGG_TCZ %>% filter(Adjusted.P.value < 0.05)

#Plot the filtered hallmark pathways
#jpeg("HallmarkPathways_TCZ_PostvsPre.jpeg", width = 2500, height = 1800, res = 300)
ggplot(Hallmark_TCZ) +
  geom_col(aes(
    x = reorder(Term, Odds.Ratio),
    y = Odds.Ratio,
    fill = -log10(Adjusted.P.value)
  ), color = "black") +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  coord_flip() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  labs(y = "Odds Ratio", 
       x = "Pathways", 
       title = "HallmarkPathway - TCZ_PostvsPre (Adj.P.Val < 0.05)",
       fill = "-log10(Adjusted.P.value)")  
#dev.off()

#Plot the filtered wiki pathways
#jpeg("WikiPathways_TCZ_PostvsPre.jpeg", width = 2500, height = 1800, res = 300)
ggplot(Wiki_TCZ %>% slice_min(Adjusted.P.value, n = 50)) +
  geom_col(aes(
    x = reorder(Term, Odds.Ratio),
    y = Odds.Ratio,
    fill = -log10(Adjusted.P.value)
  ), color = "black") +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  coord_flip() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  labs(y = "Odds Ratio", 
       x = "Pathways", 
       title = "WikiPathway - TCZ_PostvsPre (Adj.P.Val < 0.05)",
       fill = "-log10(Adjusted.P.value)")
#dev.off()

#Plot the filtered kegg pathways
#jpeg("KEGGPathways_TCZ_PostvsPre.jpeg", width = 3500, height = 2800, res = 300)
ggplot(KEGG_TCZ) +
  geom_col(aes(
    x = reorder(Term, Odds.Ratio),
    y = Odds.Ratio,
    fill = -log10(Adjusted.P.value)
  ), color = "black") +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  coord_flip() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  labs(y = "Odds Ratio", 
       x = "Pathways", 
       title = "KEGGPathway - TCZ_PostvsPre (Adj.P.Val < 0.05)",
       fill = "-log10(Adjusted.P.value)")  
#dev.off()


##-----------------------------------overlap of genes at FDR 0.05 btwn TNF vs TCZ-------------------------
TNF_genes_0.05 <- TNF_PostvsPre_FDR_0.05$gene_name # TNF genes at FDR 0.05
TCZ_genes_0.05 <- TCZ_PostvsPre_FDR_0.05$gene_name #TCZ genes at FDR 0.05

# Identify overlapping genes
overlapping_TNF_TCZ_genes <- intersect(TNF_genes_0.05, TCZ_genes_0.05)

# Number of genes in each set and overlap
length_TNF_genes_0.05 <- length(TNF_genes_0.05)
length_TCZ_genes_0.05 <- length(TCZ_genes_0.05)
length_overlap_TNF_TCZ_genes <- length(overlapping_TNF_TCZ_genes)

#jpeg("VennDiagram_TNF_TCZ_genes_0.05.jpeg", width = 3500, height = 1800, res = 300)
venn.plot <- draw.pairwise.venn(
  area1 = length_TNF_genes_0.05,
  area2 = length_TCZ_genes_0.05,
  cross.area = length_overlap_TNF_TCZ_genes,
  category = c("TNF_genes_FDR_0.05", "TCZ_genes_FDR_0.05"),
  fill = c("skyblue", "pink"),
  lty = "blank"
)
#dev.off()

# Unique genes in TNF (not in TCZ)
unique_TNF_genes <- setdiff(TNF_genes_0.05, TCZ_genes_0.05)

# Unique genes in TCZ (not in TNF)
unique_TCZ_genes <- setdiff(TCZ_genes_0.05, TNF_genes_0.05)

# Combine all unique and overlapping genes
all_genes <- unique(c(TNF_genes_0.05, TCZ_genes_0.05))

# Create a data frame for venndiagram data
venn_data <- data.frame(
  gene_name = all_genes,
  TNF = ifelse(all_genes %in% TNF_genes_0.05, 1, 0),
  TCZ = ifelse(all_genes %in% TCZ_genes_0.05, 1, 0)
)

#wb <- loadWorkbook("out_limma_study1_paired.xlsx")
#addWorksheet(wb, "Overlapped_genes_TNFvsTCZ")
#writeData(wb, "Overlapped_genes_TNFvsTCZ", venn_data)
#saveWorkbook(wb, "out_limma_study1_paired.xlsx", overwrite = T)

##-------------------------------Fisher's test on overlapped genes between TNF and TCZ treatments------------------------
#Input parameters for gene analysis
TNF_TCZ_allgenes <- length(TNF_PostvsPre$gene_name) + length(TCZ_PostvsPre$gene_name) #total no of genes analysed
setdiff_TNF_TCZ_genes <- length(TNF_genes_0.05) - length(overlapping_TNF_TCZ_genes) # Genes at FDR 0.05 in TNF but not in TCZ
setdiff_TCZ_TNF_genes <- length(TCZ_genes_0.05) - length(overlapping_TNF_TCZ_genes) # Genes at FDR 0.05 in TCZ but not in TNF
intersect_TNF_TCZ_genes <- length(overlapping_TNF_TCZ_genes)

# Calculate union and non-overlapping genes
union_TNF_TCZ_genes <- setdiff_TNF_TCZ_genes + setdiff_TCZ_TNF_genes - intersect_TNF_TCZ_genes
non_overlapping_genes <- TNF_TCZ_allgenes - union_TNF_TCZ_genes

# Create the contingency table 
contingency_table <- matrix(
  c(non_overlapping_genes, setdiff_TNF_TCZ_genes, setdiff_TCZ_TNF_genes, intersect_TNF_TCZ_genes),
  nrow = 2,
  byrow = TRUE
)

# Perform Fisher's Exact Test 
fisher_genes <- fisher.test(contingency_table, alternative = "greater")

##----------------------------------------------------------Hypergeometric test-----------------------------------------------
Hallmark <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt") # from MSigDB
Wiki <- gmtPathways("c2.cp.wikipathways.v2024.1.Hs.symbols.gmt") #from MSigDB
KEGG <- gmtPathways("c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt") #from MSigDB


TNF_DEGs_unique <- TNF_genes_0.05[!TNF_genes_0.05 %in% overlapping_TNF_TCZ_genes]
TNF_all_genes <- TNF_PostvsPre$gene_name
pathways <- c(Hallmark, Wiki, KEGG)

# Hypergeometric Test Function
hypergeom_test_TNF <- function(pathway_genes, DEGs, all_genes) {
  # Parameters for the hypergeometric test
  K <- length(pathway_genes)                  # Size of the pathway
  n <- length(DEGs)                           # Size of the DEG list
  N <- length(all_genes)                      # Total number of background genes
  x <- length(intersect(pathway_genes, DEGs)) # Overlap between DEGs and the pathway
              
  # Perform the hypergeometric test
  p_value <- phyper(q = x - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
  return(list(overlap = x, p_value = p_value))
}

# Apply the test to all pathways
TNF_results <- lapply(pathways, function(pathway_genes) {
  hypergeom_test_TNF(pathway_genes, DEGs = TNF_DEGs_unique, all_genes = TNF_all_genes)
})

# Convert results into a readable data frame
TNF_results_df <- do.call(rbind, lapply(names(TNF_results), function(pathway) {
  data.frame(
    Pathway = pathways,
    Overlap = TNF_results[[pathways]]$overlap,
    P_Value = TNF_results[[pathways]]$p_value
  )
}))

#Adjust for multiple testing using FDR
TNF_results_df$Adjusted_P_Value <- p.adjust(TNF_results_df$P_Value, method = "fdr")

# Filter significant pathways
TNF_results_df[TNF_results_df$Adjusted_P_Value < 0.05, ]
TNF_results_df[TNF_results_df$P_Value < 0.01, ] #significant pathways

## TCZ group
TCZ_DEGs_unique <- TCZ_genes_0.05[!TCZ_genes_0.05 %in% overlapping_TNF_TCZ_genes]
TCZ_all_genes <- TCZ_PostvsPre$gene_name

hypergeom_test_TCZ <- function(pathway_genes, DEGs, all_genes) {
  # Parameters for the hypergeometric test
  K <- length(pathway_genes)                  # Size of the pathway
  n <- length(DEGs)                           # Size of the DEG list
  N <- length(all_genes)                      # Total number of background genes
  x <- length(intersect(pathway_genes, DEGs)) # Overlap between DEGs and the pathway
  
  # Perform the hypergeometric test
  p_value <- phyper(q = x - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
  return(list(overlap = x, p_value = p_value))
}

# Apply the test to all pathways
TCZ_results <- lapply(pathways, function(pathway_genes) {
  hypergeom_test_TCZ(pathway_genes, DEGs = TCZ_DEGs_unique, all_genes = TCZ_all_genes)
})

TCZ_results_df <- do.call(rbind, lapply(names(TCZ_results), function(pathway) {
  data.frame(
    Pathway = pathway,
    Overlap = TCZ_results[[pathway]]$overlap,
    P_Value = TCZ_results[[pathway]]$p_value
  )
}))

#Adjust for multiple testing using FDR
TCZ_results_df$Adjusted_P_Value <- p.adjust(TCZ_results_df$P_Value, method = "fdr")

# Filter significant pathways
TCZ_results_df[TCZ_results_df$Adjusted_P_Value < 0.05, ]
TCZ_results_df[TCZ_results_df$P_Value < 0.01, ]#significant pathways
TCZ_results_df[TCZ_results_df$P_Value < 0.05, ]

##-------------------------------Enrichr- Overrepresentation analysis for unique genes in TNF and TCZ groups------------------------

TNF_unique_genes <- data.frame(gene_names = TNF_genes_0.05[!TNF_genes_0.05 %in% overlapping_TNF_TCZ_genes])
TCZ_unique_genes <- data.frame(gene_names = TCZ_genes_0.05[!TCZ_genes_0.05 %in% overlapping_TNF_TCZ_genes])

#write.xlsx(TCZ_unique_genes, "TCZ_unique_genes.xlsx")

## Over representation analysis - enrichr
Hallmark_TNF_unique <- read.delim("Hallmark_TNF_unique_genes.txt", sep = '\t')
Wiki_TNF_unique <- read.delim("Wiki_TNF_unique_genes.txt", sep = '\t')
KEGG_TNF_unique <- read.delim("KEGG_TNF_unique_genes.txt", sep = '\t')

#write.xlsx(Hallmark_TNF_unique, "out_enrichr_TNF_unique_genes.xlsx")
#wb_tnf_unque <- loadWorkbook("out_enrichr_TNF_unique_genes.xlsx")
#addWorksheet(wb_tnf_unque, "KEGGPathways")
#writeData(wb_tnf_unque, "KEGGPathways", KEGG_TNF_unique)
#saveWorkbook(wb_tnf_unque, "out_enrichr_TNF_unique_genes.xlsx", overwrite = T)

# filter the hallmark, wiki and kegg pathways based on p.val less than 0.01
Hallmark_TNF_Pval0.01 <- Hallmark_TNF_unique %>% filter(P.value < 0.01)
Wiki_TNF_Pval0.01 <- Wiki_TNF_unique %>% filter(P.value < 0.01)
KEGG_TNF_Pval0.01 <- KEGG_TNF_unique %>% filter(P.value < 0.01)

#Plot the filtered hallmark pathways
#jpeg("HallmarkPathways_TNF_PostvsPre_unique.jpeg", width = 2500, height = 1800, res = 300)
ggplot(Hallmark_TNF_Pval0.01) +
  geom_col(aes(
    x = reorder(Term, Odds.Ratio),
    y = Odds.Ratio,
    fill = -log10(P.value)
  ), color = "black") +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  coord_flip() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  labs(y = "Odds Ratio", 
       x = "Pathways", 
       title = "HallmarkPathway - TNF_PostvsPre_unique_genes
       (P.Val < 0.01)",
       fill = "-log10(P.value)")  
#dev.off()

#Plot the filtered wiki pathways
#jpeg("WikiPathways_TNF_PostvsPre_unique.jpeg", width = 2500, height = 1800, res = 300)
ggplot(Wiki_TNF_Pval0.01) +
  geom_col(aes(
    x = reorder(Term, Odds.Ratio),
    y = Odds.Ratio,
    fill = -log10(P.value)
  ), color = "black") +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  coord_flip() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  labs(y = "Odds Ratio", 
       x = "Pathways", 
       title = "WikiPathway - TNF_PostvsPre_unique_genes 
       (P.Val < 0.01)",
       fill = "-log10(P.value)")
#dev.off()

#Plot the filtered kegg pathways
#jpeg("KEGGPathways_TNF_PostvsPre_unique.jpeg", width = 3500, height = 2800, res = 300)
ggplot(KEGG_TNF_Pval0.01) +
  geom_col(aes(
    x = reorder(Term, Odds.Ratio),
    y = Odds.Ratio,
    fill = -log10(P.value)
  ), color = "black") +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  coord_flip() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  labs(y = "Odds Ratio", 
       x = "Pathways", 
       title = "KEGGPathway - TNF_PostvsPre_unique_genes 
       (P.Val < 0.01)",
       fill = "-log10(P.value)")  
#dev.off()

Hallmark_TCZ_unique <- read.delim("Hallmark_TCZ_unique_genes.txt", sep = '\t')
Wiki_TCZ_unique <- read.delim("Wiki_TCZ_unique_genes.txt", sep = '\t')
KEGG_TCZ_unique <- read.delim("KEGG_TCZ_unique_genes.txt", sep = '\t')

#write.xlsx(Hallmark_TCZ_unique, "out_enrichr_TCZ_unique_genes.xlsx")
#wb_tcz_unque <- loadWorkbook("out_enrichr_TCZ_unique_genes.xlsx")
#addWorksheet(wb_tcz_unque, "KEGGPathways")
#writeData(wb_tcz_unque, "KEGGPathways", KEGG_TCZ_unique)
#saveWorkbook(wb_tcz_unque, "out_enrichr_TCZ_unique_genes.xlsx", overwrite = T)

# filter the hallmark, wiki and kegg pathways based on p.val less than 0.01
Hallmark_TCZ_Pval0.01 <- Hallmark_TCZ_unique %>% filter(P.value < 0.01)
Wiki_TCZ_Pval0.01 <- Wiki_TCZ_unique %>% filter(P.value < 0.01)
KEGG_TCZ_Pval0.01 <- KEGG_TCZ_unique %>% filter(P.value < 0.01)

#Plot the filtered hallmark pathways
#jpeg("ORA_Hallmark_TCZ_PostvsPre_unique_DEGs.jpeg", width = 2500, height = 1800, res = 300)
ggplot(Hallmark_TCZ_Pval0.01) +
  geom_col(aes(
    x = reorder(Term, Odds.Ratio),
    y = Odds.Ratio,
    fill = -log10(P.value)
  ), color = "black") +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  coord_flip() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  labs(y = "Odds Ratio", 
       x = "Pathways", 
       title = "HallmarkPathway - TCZ_PostvsPre_unique_DEGs
       (Pval < 0.01)",
       fill = "-log10(P.value)")  
#dev.off()

#Plot the filtered wiki pathways
#jpeg("ORA_Wiki_TCZ_PostvsPre_unique_DEGs.jpeg", width = 2500, height = 1800, res = 300)
ggplot(Wiki_TCZ_Pval0.01) +
  geom_col(aes(
    x = reorder(Term, Odds.Ratio),
    y = Odds.Ratio,
    fill = -log10(P.value)
  ), color = "black") +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  coord_flip() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  labs(y = "Odds Ratio", 
       x = "Pathways", 
       title = "WikiPathway - TCZ_PostvsPre_unique_DEGs
       (Pval < 0.01)",
       fill = "-log10(P.value)")
#dev.off()

#Plot the filtered kegg pathways
#jpeg("ORA_KEGG_TCZ_PostvsPre_unique_DEGs.jpeg", width = 3500, height = 2800, res = 300)
ggplot(KEGG_TCZ_Pval0.01) +
  geom_col(aes(
    x = reorder(Term, Odds.Ratio),
    y = Odds.Ratio,
    fill = -log10(P.value)
  ), color = "black") +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  coord_flip() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  labs(y = "Odds Ratio", 
       x = "Pathways", 
       title = "KEGGPathway - TCZ_PostvsPre_unique_DEGs 
       (Pval < 0.01)",
       fill = "-log10(P.value)")  
#dev.off()



##-----------------------------overlap of pathways found at Pval 0.01 between TNF vs TCZ groups -----------------------
TNF_Hallmark_pathways <- unique(Hallmark_TNF_Pval0.01$Term)
TNF_Wiki_pathways <- unique(Wiki_TNF_Pval0.01$Term)
TNF_KEGG_pathways <- unique(KEGG_TNF_Pval0.01$Term)

TCZ_Hallmark_pathways <- unique(Hallmark_TCZ_Pval0.01$Term)
TCZ_Wiki_pathways <- unique(Wiki_TCZ_Pval0.01$Term)
TCZ_KEGG_pathways <- unique(KEGG_TCZ_Pval0.01$Term)

#store hallmark pathway lists in a named list
hl_pathways <- list(
  Hallmark_TNF = TNF_Hallmark_pathways,
  Hallmark_TCZ = TCZ_Hallmark_pathways
)

# visualize overlapping pathways at FDR 0.05
venn.plot.hl <- venn.diagram(
  x = hl_pathways,
  category.names = names(hl_pathways),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.3,
  cat.cex = 1.2,
  cat.fontface = "bold",
  cex = 1.5,
  scaled= F,
  margin = 0.1,
  height = 2000,
  width = 2000,
  cat.pos = c(0,0),     
  cat.dist = c(0.1, 0.1)
)

#jpeg("VennDiagram_Hallmark_unique_DEGs.jpeg", width = 2500, height = 1800, res = 300 )
grid.draw(venn.plot.hl)
#dev.off()

TNF_Hallmark_pathways <- data.frame(Term = TNF_Hallmark_pathways)
TNF_Hallmark_pathways <- TNF_Hallmark_pathways %>% select(Term) %>% mutate(Hallmark_TNF = 1)
TCZ_Hallmark_pathways <- data.frame(Term = TCZ_Hallmark_pathways)
TCZ_Hallmark_pathways <- TCZ_Hallmark_pathways %>% select(Term) %>% mutate(Hallmark_TCZ = 1)

combine_hl <- full_join(TNF_Hallmark_pathways, TCZ_Hallmark_pathways, by = "Term") %>%
  replace(is.na(.), 0)

#write.xlsx(combine_hl, "Overlapping_enrichr_pathways_uniquegenes.xlsx")

 #store wiki pathway lists in a named list
wk_pathways <- list(
  Wiki_TNF = TNF_Wiki_pathways,
  Wiki_TCZ = TCZ_Wiki_pathways
)

#  visualize overlapping pathways at FDR 0.05
venn.plot.wk <- venn.diagram(
  x = wk_pathways,
  category.names = names(wk_pathways),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.3,
  cat.cex = 1.2,
  cat.fontface = "bold",
  cex = 1.5,
  scaled= F,
  margin = 0.1,
  height = 2000,
  width = 2000,
  cat.pos = c(0,0),     
  cat.dist = c(0.1, 0.1)
)

#jpeg("VennDiagram_Wiki_unique_DEGs.jpeg", width = 2500, height = 1800, res = 300 )
grid.draw(venn.plot.wk)
#dev.off()

TNF_Wiki_pathways <- data.frame(Term = TNF_Wiki_pathways)
TNF_Wiki_pathways <- TNF_Wiki_pathways %>% select(Term) %>% mutate(Wiki_TNF = 1)
TCZ_Wiki_pathways <- data.frame(Term = TCZ_Wiki_pathways)
TCZ_Wiki_pathways <- TCZ_Wiki_pathways %>% select(Term) %>% mutate(Wiki_TCZ = 1)

combine_wk <- full_join(TNF_Wiki_pathways, TCZ_Wiki_pathways, by = "Term") %>%
  replace(is.na(.), 0)

#store kegg pathway lists in a named list
Kg_pathways <- list(
  KEGG_TNF = TNF_KEGG_pathways,
  KEGG_TCZ = TCZ_KEGG_pathways
)

#  visualize overlapping pathways at FDR 0.05
venn.plot.Kg <- venn.diagram(
  x = Kg_pathways,
  category.names = names(Kg_pathways),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.3,
  cat.cex = 1.2,
  cat.fontface = "bold",
  cex = 1.5,
  scaled= F,
  margin = 0.1,
  height = 2000,
  width = 2000,
  cat.pos = c(0,0),     
  cat.dist = c(0.1, 0.1)
)

#jpeg("VennDiagram_KEGG_unique_DEGs.jpeg", width = 2500, height = 1800, res = 300 )
grid.draw(venn.plot.Kg)
#dev.off()

TNF_KEGG_pathways <- data.frame(Term = TNF_KEGG_pathways)
TNF_KEGG_pathways <- TNF_KEGG_pathways %>% select(Term) %>% mutate(KEGG_TNF = 1)
TCZ_KEGG_pathways <- data.frame(Term = TCZ_KEGG_pathways)
TCZ_KEGG_pathways <- TCZ_KEGG_pathways %>% select(Term) %>% mutate(KEGG_TCZ = 1)

combine_Kg <- full_join(TNF_KEGG_pathways, TCZ_KEGG_pathways, by = "Term") %>%
  replace(is.na(.), 0)

#wb_pathwys <- loadWorkbook("Overlapping_enrichr_pathways_uniquegenes.xlsx")
#addWorksheet(wb_pathwys, "KEGGPathways")
#writeData(wb_pathwys, "KEGGPathways", combine_Kg)
#saveWorkbook(wb_pathwys, "overlapping_enrichr_pathways_uniquegenes.xlsx", overwrite = T)




##-----------------------------------------------------study2-------------------------------------------------------------------------------
#check if there are any duplicated genes in raw count data 
any(duplicated(study2_raw_counts$gene_name)) #found 1624 duplicated genes

S2_dup_genes <- study2_raw_counts[duplicated(study2_raw_counts$gene_name), ] # add the duplicated genes data into a new df
S2_dup_genes <-  S2_dup_genes %>%                          # concatenate the duplicated gene name with its corresponding ensembl_id
  mutate(gene_name = paste0(ensembl_gene_id, "_", gene_name))

#replace the modified duplicated gene name values in the raw counts data
study2_raw_counts$gene_name <- ifelse(
  rownames(study2_raw_counts) %in% rownames(S2_dup_genes),  # Check if row names match
  S2_dup_genes$gene_name[match(rownames(study2_raw_counts), rownames(S2_dup_genes))],  # Replace with modified names
  study2_raw_counts$gene_name  # Keep original name if no match
)

## create a DGE list for study2 data
study2_raw_counts <- study2_raw_counts[,-1]
study2_raw_counts <- subset(study2_raw_counts, select = -c(WD41,WD42, WD43))
study2_DGE <- DGEList(counts = study2_raw_counts)

## TMM normalization
study2_DGE <- calcNormFactors(study2_DGE)

## log2 transformation of count data and prior count = 3
study2_logCPM <- cpm(study2_DGE$counts, log = T, prior.count = 3) # logcpm values without gene names

study2_logCPM_file <- data.frame(Gene = study2_DGE$genes$gene_name, study2_logCPM) #logcpm values with gene names


## add group data to DGE list
study2_groups <- study2_sample_data %>% select(c(1,5))
study2_groups <- study2_groups[-c(19,20,21),]
rownames(study2_groups) <- NULL
study2_groups$Excerise_type <- ifelse(
  grepl("EX", study2_groups$SAMPLE.ID), "EX",
  ifelse(
    grepl("AFL", study2_groups$SAMPLE.ID), "AFL",
    ifelse(
      grepl("SED", study2_groups$SAMPLE.ID), "SED", 
      NA)))

study2_groups$Patient <- ifelse(
  grepl("ABC", study2_groups$SAMPLE.ID), "Sub1",
  ifelse(
   grepl("MAS", study2_groups$SAMPLE.ID), "Sub2",
    ifelse(
      grepl("HCG", study2_groups$SAMPLE.ID), "Sub3",
      ifelse(
        grepl("RSP", study2_groups$SAMPLE.ID), "Sub4",
        ifelse(
          grepl("RSA", study2_groups$SAMPLE.ID), "Sub5",
          ifelse(
            grepl("MKA", study2_groups$SAMPLE.ID), "Sub6",
              ifelse(
                grepl("ML", study2_groups$SAMPLE.ID), "Sub7",
                ifelse(
                  grepl("RBM", study2_groups$SAMPLE.ID), "Sub8",
                  ifelse(
                    grepl("MSS", study2_groups$SAMPLE.ID), "Sub9",
                    NA)))))))))

study2_groups$Patient_Excerise <-  paste(study2_groups$Patient,study2_groups$Raw.count.file, sep = '_') #combine patient and excerise data together

## add the group values to DGE list
study2_DGE$samples$group <- study2_groups$Excerise_type

# Filter the lowly expressed genes
S2_keep <- rowSums(cpm(study2_DGE) >1) >= 13
study2_DGE <- study2_DGE[S2_keep,,keep.lib.sizes = F]
length(rownames(study2_DGE$counts)) #13112 genes remained after filtering

## Density plot
S2_L <- mean(study2_DGE$samples$lib.size) * 1e-6
S2_M <- median(study2_DGE$samples$lib.size) * 1e-6
S2_logCPM_thres <- log2(10/S2_M + 2/S2_L)

#jpeg("Study2_DensityPlot.jpeg", width = 2500, height = 1800, res = 300)
#Density Plot 
S2_samples <- ncol(study2_DGE)
col <- brewer.pal(S2_samples, "Dark2")
par(mfrow=c(1,2))
plot(density(study2_logCPM[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main=" Raw_data", xlab="Log-cpm")
abline(v= S2_logCPM_thres, lty=3)
for (i in 2:S2_samples){
  den <- density(study2_logCPM[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(study2_DGE), text.col=col, bty="n")
study2_logCPM <- cpm(study2_DGE, log=TRUE, prior.count = 3)
plot(density(study2_logCPM[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main=" Filtered_data", xlab="Log-cpm")
abline(v=S2_logCPM_thres, lty=3)
for (i in 2:S2_samples){
  den <- density(study2_logCPM[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(study2_DGE), text.col=col, bty="n")
#dev.off()


## create MDS plot to check the variability in gene expr data
S2_group_colors <- ifelse(
  study2_DGE$samples$group == "SED", "blue",
  ifelse(
    study2_DGE$samples$group == "AFL", "darkgreen",
    ifelse(study2_DGE$samples$group == "EX", "red", "black")  # Default to black for any unmatched group
  )
)

#jpeg("Study2_MDSplot.jpeg", width = 2500, height = 1800, res = 300)
plotMDS(study2_logCPM, labels = study2_DGE$samples$group, col = S2_group_colors, cex = 1)
#dev.off()

## PCA plot 
## Principal Component Analysis 
S2_data <- data.table(study2_logCPM) # convert logcpm matrix to a data table obj
S2_tdata <- data.table::transpose(S2_data) # transpose the data table
colnames(S2_tdata) <- study2_DGE$genes$gene_name #set colnames of the transposed data
rownames(S2_tdata) <- colnames(study2_logCPM) #set rownames of the transposed data

S2_pca <- prcomp(S2_tdata, scale. = T) # pca on transposed data

##screeplot 
screeplot(S2_pca)

#dev.off()

## proportion of variance explained by each component
summary(S2_pca)

S2_pca_data <- data.frame(S2_pca$x) # convert principal component values to a dataframe obj
S2_pca_data$Exercise_type <- study2_DGE$samples$group #Add the groups column to pca dataframe
S2_pca_data$Patient <- study2_groups$Patient_Excerise

# Assign colors to groups
S2_pca_groups <- unique(study2_DGE$samples$group)
S2_colors <- brewer.pal(length(S2_pca_groups), "Set2")
S2_pca_data$Color <- S2_colors[as.numeric(factor(S2_pca_data$Exercise_type, levels = S2_pca_groups))]

#jpeg("Study2_PCA.jpeg", width = 2500, height = 1800, res = 300)
# Plot PCA
ggplot(data= S2_pca_data, aes(x=PC1, y=PC2)) +
  geom_point(size=5, aes(color=Exercise_type, shape = Exercise_type))+
  geom_hline(yintercept=0, colour="black", linetype="dashed") +
  geom_vline(xintercept=0, colour="black", linetype="dashed") +
  geom_text_repel(label= S2_pca_data$Patient, colour="black", size=3, nudge_x=0.1, nudge_y=0.1, max.overlaps = Inf) +  
  xlab("PC1 (16%)")+
  ylab("PC2(13%)")+ 
  theme(axis.text=element_text(size=14, colour="black"),axis.title=element_text(size=14, face="bold", colour="black")) + 
  theme(panel.background = element_rect(fill = "white", colour = 'black'))+
  ggtitle("Study2_PCA_plot_Wagner_Study")
#dev.off()

##Differential gene expression analysis - Paired data
study2_groups <- study2_groups %>% group_by(Patient) %>% mutate(Patient_Combined = paste0(Patient, "_", paste(unique(Raw.count.file), collapse = "_"))) %>%
                     ungroup()
S2_subject <- factor(study2_groups$Patient_Combined)
S2_subject <- factor(S2_subject, levels = unique(S2_subject)) #get the subject data#

S2_exercise_type <- factor(study2_DGE$samples$group)
S2_exercise_type <- factor(S2_exercise_type, levels = unique(S2_exercise_type))

S2_paired_data <- data.frame(Samplenames = colnames(study2_DGE$counts), Exercise_type = S2_exercise_type, Subject = S2_subject)

## create a design matrix 
S2_design <- model.matrix(~0+Exercise_type+Subject, data = S2_paired_data)
colnames(S2_design) <- gsub("Exercise_type", "", colnames(S2_design))
colnames(S2_design) <- gsub("Subject", "", colnames(S2_design))
rownames(S2_design) <- S2_paired_data$Samplenames

## contrast matrix 
S2_contr.mtrx <- makeContrasts(
  EX_vs_SED = EX - SED,
  EX_vs_AFL = EX -AFL,
  AFL_vs_SED = AFL - SED,
  levels = S2_design)
S2_contr.mtrx

## run voom
S2_v <- voom(study2_DGE, S2_design, plot = T)
S2_fit <- lmFit(S2_v, S2_design)
S2_fit <- contrasts.fit(S2_fit, S2_contr.mtrx)
S2_efit <- eBayes(S2_fit)

EX_vs_SED <- topTable(S2_efit, coef = 1, n = Inf, adjust.method = "fdr", sort.by = "P")
EX_vs_AFL <- topTable(S2_efit, coef = 2, n = Inf, adjust.method = "fdr", sort.by = "P")
AFL_vs_SED <- topTable(S2_efit, coef = 3, n = Inf, adjust.method = "fdr", sort.by = "P")

EX_vs_SED_Pval_0.01 <- EX_vs_SED %>% filter(P.Value < 0.01)
EX_vs_AFL_Pval_0.01 <- EX_vs_AFL %>% filter(P.Value < 0.01)
AFL_vs_SED_Pval_0.01 <- AFL_vs_SED %>% filter(P.Value < 0.01)



##--------------------------------------------------------Removal of Outliers-----------------------------------
study2_raw_counts <- subset(study2_raw_counts, select = -c(WD38,WD39, WD40)) #remove the outlier data from the raw counts
study2_DGE_out <- DGEList(counts = study2_raw_counts)

## TMM normalization
study2_DGE_out <- calcNormFactors(study2_DGE_out)

## log2 transformation of count data and prior count = 3
study2_logCPM_out <- cpm(study2_DGE_out$counts, log = T, prior.count = 3) # logcpm values without gene names

study2_logCPM_out_file <- data.frame(Gene = study2_DGE_out$genes$gene_name, study2_logCPM_out) #logcpm values with gene names
#write.xlsx(study2_logCPM_out_file, "study2_logCPM.xlsx")

study2_groups_out <- study2_groups[-c(16,17,18),] #remove the outlier data from the sample information
study2_groups_out <- study2_groups_out %>%
  mutate(Patient = paste0("Sub", ceiling(row_number() / 3))) #adjust the patient data based on the outliers removed
study2_groups_out$Patient_Excerise <- paste0("Sub", ceiling(seq_len(nrow(study2_groups_out)) /3), "_", sub(".*_(WD.*)", "\\1", study2_groups_out$Patient_Excerise)) #same with patient_exercise data
study2_groups_out$Patient_Combined <- paste0("Sub", ceiling(seq_len(nrow(study2_groups_out)) / 3),  "_", sub(".*?_(WD.*)", "\\1", study2_groups_out$Patient_Combined)) #same with the  paired data

# add the new group values to DGE list
study2_DGE_out$samples$group <- study2_groups_out$Excerise_type

# Filter the lowly expressed genes
S2_keep1 <- rowSums(cpm(study2_DGE_out) >1) >= 12
study2_DGE_out <- study2_DGE_out[S2_keep1,,keep.lib.sizes = F]
length(rownames(study2_DGE_out$counts)) #12989 genes remained after filtering

## Density plot
S2_L1 <- mean(study2_DGE_out$samples$lib.size) * 1e-6
S2_M1 <- median(study2_DGE_out$samples$lib.size) * 1e-6
S2_logCPM_thres1 <- log2(10/S2_M1 + 2/S2_L1)

#jpeg("Study2_DensityPlot_without_outliers.jpeg", width = 2500, height = 1800, res = 300)
#Density Plot 
S2_samples1 <- ncol(study2_DGE_out)
col <- brewer.pal(S2_samples1, "Dark2")
par(mfrow=c(1,2))
plot(density(study2_logCPM_out[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main=" Raw_data", xlab="Log-cpm")
abline(v= S2_logCPM_thres1, lty=3)
for (i in 2:S2_samples1){
  den <- density(study2_logCPM_out[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(study2_DGE_out), text.col=col, bty="n")
study2_logCPM_out <- cpm(study2_DGE_out, log=TRUE, prior.count = 3)
plot(density(study2_logCPM_out[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main=" Filtered_data", xlab="Log-cpm")
abline(v=S2_logCPM_thres1, lty=3)
for (i in 2:S2_samples1){
  den <- density(study2_logCPM_out[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(study2_DGE_out), text.col=col, bty="n")
#dev.off()

## PCA plot 
## Principal Component Analysis 
S2_data1 <- data.table(study2_logCPM_out) # convert logcpm matrix to a data table obj
S2_tdata1 <- data.table::transpose(S2_data1) # transpose the data table
colnames(S2_tdata1) <- study2_DGE_out$genes$gene_name #set colnames of the transposed data
rownames(S2_tdata1) <- colnames(study2_logCPM_out) #set rownames of the transposed data

S2_pca_out <- prcomp(S2_tdata1, scale. = T) # pca on transposed data

##screeplot 
screeplot(S2_pca_out)

#dev.off()

## proportion of variance explained by each component
summary(S2_pca_out)

S2_pca_data_out <- data.frame(S2_pca_out$x) # convert principal component values to a dataframe obj
S2_pca_data_out$Exercise_type <- study2_DGE_out$samples$group #Add the groups column to pca dataframe
S2_pca_data_out$Patient <- study2_groups_out$Patient_Excerise

# Assign colors to groups
S2_pca_groups_out <- unique(study2_DGE_out$samples$group)
S2_colors_out <- brewer.pal(length(S2_pca_groups_out), "Set2")
S2_pca_data_out$Color <- S2_colors_out[as.numeric(factor(S2_pca_data_out$Exercise_type, levels = S2_pca_groups_out))]

#jpeg("Study2_PCA_without_outliers.jpeg", width = 2500, height = 1800, res = 300)
# Plot PCA
ggplot(data= S2_pca_data_out, aes(x=PC1, y=PC2)) +
  geom_point(size=5, aes(color=Exercise_type, shape = Exercise_type))+
  geom_hline(yintercept=0, colour="black", linetype="dashed") +
  geom_vline(xintercept=0, colour="black", linetype="dashed") +
  geom_text_repel(label= S2_pca_data_out$Patient, colour="black", size=3, nudge_x=0.1, nudge_y=0.1, max.overlaps = Inf) +  
  xlab("PC1 (15%)")+
  ylab("PC2(12%)")+ 
  theme(axis.text=element_text(size=14, colour="black"),axis.title=element_text(size=14, face="bold", colour="black")) + 
  theme(panel.background = element_rect(fill = "white", colour = 'black'))+
  ggtitle("Study2_PCA_plot_Wagner_Study")
#dev.off()

##Differential expression analysis - after outlier removal
S2_subject_out <- factor(study2_groups_out$Patient_Combined)
S2_subject_out <- factor(S2_subject_out, levels = unique(S2_subject_out)) #get the subject data#

S2_exercise_type_out <- factor(study2_DGE_out$samples$group)
S2_exercise_type_out <- factor(S2_exercise_type_out, levels = unique(S2_exercise_type_out))

S2_paired_data_out <- data.frame(Samplenames = colnames(study2_DGE_out$counts), Exercise_type = S2_exercise_type_out, Subject = S2_subject_out)

## create a design matrix 
S2_design_out <- model.matrix(~0+Exercise_type+Subject, data = S2_paired_data_out)
colnames(S2_design_out) <- gsub("Exercise_type", "", colnames(S2_design_out))
colnames(S2_design_out) <- gsub("Subject", "", colnames(S2_design_out))
rownames(S2_design_out) <- S2_paired_data_out$Samplenames

## contrast matrix 
S2_contr.mtrx_out <- makeContrasts(
  EX_vs_SED = EX - SED,
  EX_vs_AFL = EX -AFL,
  AFL_vs_SED = AFL - SED,
  levels = S2_design_out)
S2_contr.mtrx_out

#run trend
trendfit <- lmFit(study2_logCPM_out, S2_design_out)
trendfit <- contrasts.fit(trendfit, S2_contr.mtrx_out)
efit <- eBayes(trendfit, trend = T)

EX_vs_SED_out <- topTable(efit, coef = 1, n = Inf, adjust.method = "fdr", sort.by = "P")
EX_vs_AFL_out <- topTable(efit, coef = 2, n = Inf, adjust.method = "fdr", sort.by = "P")
AFL_vs_SED_out <- topTable(efit, coef = 3, n = Inf, adjust.method = "fdr", sort.by = "P")

# add gene names 
EX_vs_SED_out <- EX_vs_SED_out %>% mutate(gene_name = study2_logCPM_out_file[rownames(EX_vs_SED_out), "Gene"]) %>% select(gene_name, everything())
EX_vs_AFL_out <- EX_vs_AFL_out %>% mutate(gene_name = study2_logCPM_out_file[rownames(EX_vs_AFL_out), "Gene"]) %>% select(gene_name, everything())
AFL_vs_SED_out <- AFL_vs_SED_out %>% mutate(gene_name = study2_logCPM_out_file[rownames(AFL_vs_SED_out), "Gene"]) %>% select(gene_name, everything())

#write.xlsx(EX_vs_SED_out, "out_limma_exvssed.xlsx")
#write.xlsx(EX_vs_AFL_out, "out_limma_exvsafl.xlsx")
#write.xlsx(AFL_vs_SED_out, "out_limma_aflvssed.xlsx")

#volcano plot for Ex_vs_SED
#jpeg("Volcanoplot_EX_vs_SED.jpeg", width = 2500, height = 1800, res = 300)
ggplot(EX_vs_SED_out, aes(x= logFC ,
                          y=-log10(P.Value), 
                          color= -log10(P.Value)))+
  geom_point(size=2)+
  geom_point(data=subset(EX_vs_SED_out %>% slice_min(P.Value, n=20)),
             aes(x= logFC, 
                 y=-log10(P.Value), 
                 color= -log10(P.Value)))+
  scale_color_viridis(option="rocket", direction=1)+
  geom_text_repel(data=subset(EX_vs_SED_out %>% slice_min(P.Value, n=20)),
                  aes(x=logFC, 
                      y=-log10(P.Value), 
                      color= -log10(P.Value),label = gene_name), 
                  color = "black", size=4, box.padding = unit(0.2, "lines"),segment.color="gray25",
                  segment.size=0.25, point.padding = unit(0.2, "lines"), max.overlaps = Inf)+
  geom_hline(yintercept=-log10(0.01), linetype="dashed") +
  geom_vline(xintercept=c(1,-1),linetype="dashed")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text (size = 10, face ="bold", colour = "black"))+
  theme(axis.title = element_text (size = 18, face = "bold", colour = "black")) +
  theme(strip.text = element_text(size=10, face="bold", color = "black"))+
  ggtitle("volcano_Plot_Top20genes_EXvsAFL")+
  theme(plot.title=element_text(face="bold"))+
  xlab("log2FC (EX vs. SED)") +
  ylab("-log10 P.Value (EX vs. SED)")+
  theme(legend.position="bottom") #to maximize graph space by removing legends
#dev.off()

#jpeg("Volcanoplot_EX_vs_AFL.jpeg", width = 2500, height = 1800, res = 300)
ggplot(EX_vs_AFL_out, aes(x= logFC ,
                          y=-log10(P.Value), 
                          color= -log10(P.Value)))+
  geom_point(size=2)+
  geom_point(data=subset(EX_vs_AFL_out %>% slice_min(P.Value, n=20)),
             aes(x= logFC, 
                 y=-log10(P.Value), 
                 color= -log10(P.Value)))+
  scale_color_viridis(option="rocket", direction=1)+
  geom_text_repel(data=subset(EX_vs_AFL_out %>% slice_min(P.Value, n=20)),
                  aes(x=logFC, 
                      y=-log10(P.Value), 
                      color= -log10(P.Value),label = gene_name), 
                  color = "black", size=4, box.padding = unit(0.2, "lines"),segment.color="gray25",
                  segment.size=0.25, point.padding = unit(0.2, "lines"), max.overlaps = Inf)+
  geom_hline(yintercept=-log10(0.01), linetype="dashed") +
  geom_vline(xintercept=c(1,-1),linetype="dashed")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text (size = 10, face ="bold", colour = "black"))+
  theme(axis.title = element_text (size = 18, face = "bold", colour = "black")) +
  theme(strip.text = element_text(size=10, face="bold", color = "black"))+
  ggtitle("volcano_Plot_Top20genes_EXvsAFL")+
  theme(plot.title=element_text(face="bold"))+
  xlab("log2FC (EX vs. AFL)") +
  ylab("-log10 P.Value (EX vs. AFL)")+
  theme(legend.position="bottom") #to maximize graph space by removing legends
#dev.off()

#jpeg("Volcanoplot_AFL_vs_SED.jpeg", width = 2500, height = 1800, res = 300)
ggplot(AFL_vs_SED_out, aes(x= logFC ,
                          y=-log10(P.Value), 
                          color= -log10(P.Value)))+
  geom_point(size=2)+
  geom_point(data=subset(AFL_vs_SED_out %>% slice_min(P.Value, n=20)),
             aes(x= logFC, 
                 y=-log10(P.Value), 
                 color= -log10(P.Value)))+
  scale_color_viridis(option="rocket", direction=1)+
  geom_text_repel(data=subset(AFL_vs_SED_out %>% slice_min(P.Value, n=20)),
                  aes(x=logFC, 
                      y=-log10(P.Value), 
                      color= -log10(P.Value),label = gene_name), 
                  color = "black", size=4, box.padding = unit(0.2, "lines"),segment.color="gray25",
                  segment.size=0.25, point.padding = unit(0.2, "lines"), max.overlaps = Inf)+
  geom_hline(yintercept=-log10(0.01), linetype="dashed") +
  geom_vline(xintercept=c(1,-1),linetype="dashed")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text (size = 10, face ="bold", colour = "black"))+
  theme(axis.title = element_text (size = 18, face = "bold", colour = "black")) +
  theme(strip.text = element_text(size=10, face="bold", color = "black"))+
  ggtitle("volcano_Plot_Top20genes_AFLvsSED")+
  theme(plot.title=element_text(face="bold"))+
  xlab("log2FC (AFL vs. SED)") +
  ylab("-log10 P.Value (AFL vs. SED)")+
  theme(legend.position="bottom") #to maximize graph space by removing legends
#dev.off()

##Heatmaps for Ex_vs_SED comparison
EXvsSED_filtered <- EX_vs_SED_out[EX_vs_SED_out$P.Value < 0.01 & abs(EX_vs_SED_out$logFC) > 1, ]
EXvsSED.i.1 <- study2_DGE_out$genes$gene_name[which(study2_DGE_out$genes$gene_name %in% EXvsSED_filtered$gene_name)] 

EXvsSED_logCPM <- study2_logCPM_out_file[, c(TRUE, colnames(study2_logCPM_out_file)[-1] %in% study2_groups_out$Raw.count.file[study2_groups_out$Excerise_type %in% c("SED", "EX")])]

data_ExvsSED <- EXvsSED_logCPM[EXvsSED_logCPM$Gene %in% EXvsSED.i.1, -1]
rownames(data_ExvsSED) <- EXvsSED.i.1

ExvsSED.samples.group <- paste(colnames(data_ExvsSED), study2_groups_out$Excerise_type[study2_groups_out$Excerise_type %in% c("SED", "EX")], sep = "_")
colnames(data_ExvsSED) <- ExvsSED.samples.group

SED_samples <- grep("SED", colnames(data_ExvsSED), value = TRUE)  # Extract sedentary samples
EX_samples <- grep("EX", colnames(data_ExvsSED), value = T) #extract high intensity workout samples

# Reorder the columns such that sedentary and high intensity workout samples are grouped
EXvsSED_ordered_samples <- c(SED_samples, EX_samples)

#color palette
palette <- colorRampPalette(c("blue", "white", "red"))(100)

data_ExvsSED <- data_ExvsSED[, EXvsSED_ordered_samples]

# create an annotation dataframe for the samples
annot_col_ExvsSED <- data.frame(workouts = substr(colnames(data_ExvsSED),6,8))
rownames(annot_col_ExvsSED) <- colnames(data_ExvsSED)


annot_colors_EXvsSED <- list(
  workouts = c(SED = "darkgreen", EX = "darkred")
)


#jpeg("Heatmap_ExvsSED.jpeg", width = 2500, height = 1800, res = 300)
#generate heatmap
pheatmap(data_ExvsSED, show_rownames = T, scale = "row", cluster_cols= F, cluster_rows = T, color = palette, border_color = NA,
         fontsize_row = 8 , fontsize_col = 10,annotation_col = annot_col_ExvsSED, annotation_colors = annot_colors_EXvsSED, 
         annotation_legend = T, main = "Topgenes_EXvsSED_PVal_0.01")
#dev.off()


##Heatmaps for Ex_vs_AFL comparison
EXvsAFL_filtered <- EX_vs_AFL_out[EX_vs_AFL_out$P.Value < 0.01 & abs(EX_vs_AFL_out$logFC) > 1, ]
EXvsAFL.i.1 <- study2_DGE_out$genes$gene_name[which(study2_DGE_out$genes$gene_name %in% EXvsAFL_filtered$gene_name)] 

EXvsAFL_logCPM <- study2_logCPM_out_file[, c(TRUE, colnames(study2_logCPM_out_file)[-1] %in% study2_groups_out$Raw.count.file[study2_groups_out$Excerise_type %in% c("AFL", "EX")])]

data_EXvsAFL <- EXvsAFL_logCPM[EXvsAFL_logCPM$Gene %in% EXvsAFL.i.1, -1]
rownames(data_EXvsAFL) <- EXvsAFL.i.1

EXvsAFL.samples.group <- paste(colnames(data_EXvsAFL), study2_groups_out$Excerise_type[study2_groups_out$Excerise_type %in% c("AFL", "EX")], sep = "_")
colnames(data_EXvsAFL) <- EXvsAFL.samples.group

AFL_samples <- grep("AFL", colnames(data_EXvsAFL), value = TRUE)  # Extract moderate workout samples
EX_samples <- grep("EX", colnames(data_EXvsAFL), value = T) #extract high intensity workout samples

# Reorder the columns such that sedentary and high intensity workout samples are grouped
EXvsAFL_ordered_samples <- c(AFL_samples, EX_samples)

#color palette
palette <- colorRampPalette(c("blue", "white", "red"))(100)

data_EXvsAFL <- data_EXvsAFL[, EXvsAFL_ordered_samples]

# create an annotation dataframe for the samples
annot_col_EXvsAFL <- data.frame(workouts = substr(colnames(data_EXvsAFL),6,8))
rownames(annot_col_EXvsAFL) <- colnames(data_EXvsAFL)


annot_colors_EXvsAFL <- list(
  workouts = c(AFL = "darkgreen", EX = "darkred")
)


#jpeg("Heatmap_EXvsAFL.jpeg", width = 2500, height = 1800, res = 300)
#generate heatmap
pheatmap(data_EXvsAFL, show_rownames = T, scale = "row", cluster_cols= F, cluster_rows = T, color = palette, border_color = NA,
         fontsize_row = 8 , fontsize_col = 10,annotation_col = annot_col_EXvsAFL, annotation_colors = annot_colors_EXvsAFL, 
         annotation_legend = T, main = "Topgenes_EXvsAFL_PVal_0.01")
#dev.off()

##Heatmaps for AFL_vs_SED comparison
AFLvsSED_filtered <- AFL_vs_SED_out[AFL_vs_SED_out$P.Value < 0.01 & abs(AFL_vs_SED_out$logFC) > 1, ]
AFLvsSED.i.1 <- study2_DGE_out$genes$gene_name[which(study2_DGE_out$genes$gene_name %in% AFLvsSED_filtered$gene_name)] 

AFLvsSED_logCPM <- study2_logCPM_out_file[, c(TRUE, colnames(study2_logCPM_out_file)[-1] %in% study2_groups_out$Raw.count.file[study2_groups_out$Excerise_type %in% c("SED", "AFL")])]

data_AFLvsSED <- AFLvsSED_logCPM[AFLvsSED_logCPM$Gene %in% AFLvsSED.i.1, -1]
rownames(data_AFLvsSED) <- AFLvsSED.i.1

AFLvsSED.samples.group <- paste(colnames(data_AFLvsSED), study2_groups_out$Excerise_type[study2_groups_out$Excerise_type %in% c("SED", "AFL")], sep = "_")
colnames(data_AFLvsSED) <- AFLvsSED.samples.group

SED_samples <- grep("SED", colnames(data_AFLvsSED), value = TRUE)  # Extract sedentary workout samples
AFL_samples <- grep("AFL", colnames(data_AFLvsSED), value = T) #extract moderate workout samples

# Reorder the columns such that sedentary and high intensity workout samples are grouped
AFLvsSED_ordered_samples <- c(SED_samples, AFL_samples)

#color palette
palette <- colorRampPalette(c("blue", "white", "red"))(100)

data_AFLvsSED <- data_AFLvsSED[, AFLvsSED_ordered_samples]

# create an annotation dataframe for the samples
annot_col_AFLvsSED <- data.frame(workouts = substr(colnames(data_AFLvsSED),6,8))
rownames(annot_col_AFLvsSED) <- colnames(data_AFLvsSED)


annot_colors_AFLvsSED <- list(
  workouts = c(SED = "darkgreen", AFL = "darkred")
)


#jpeg("Heatmap_AFLvsSED.jpeg", width = 2500, height = 1800, res = 300)
#generate heatmap
pheatmap(data_AFLvsSED, show_rownames = T, scale = "row", cluster_cols= F, cluster_rows = T, color = palette, border_color = NA,
         fontsize_row = 8 , fontsize_col = 10,annotation_col = annot_col_AFLvsSED, annotation_colors = annot_colors_AFLvsSED, 
         annotation_legend = T, main = "Topgenes_AFLvsSED_PVal_0.01")
#dev.off()


##-------------------------------------------------New analysis - 2/10/2025 - Venndiagram of study2 genes------------------
EX_vs_SED_out_genes_0.01 <- EX_vs_SED_out %>% filter(P.Value < 0.01) %>% pull(gene_name)
EX_vs_AFL_out_genes_0.01 <- EX_vs_AFL_out %>% filter(P.Value < 0.01) %>% pull(gene_name)
AFL_vs_SED_out_genes_0.01 <- AFL_vs_SED_out %>% filter(P.Value < 0.01) %>% pull(gene_name)
overlapping_study2_genes <- purrr::reduce(list(EX_vs_SED_out_genes_0.01, EX_vs_AFL_out_genes_0.01, AFL_vs_SED_out_genes_0.01), intersect)


length_EX_vs_SED_out_0.01 <- length(EX_vs_SED_out_genes_0.01)
length_EX_vs_AFL_out_0.01 <- length(EX_vs_AFL_out_genes_0.01)
length_AFL_vs_SED_out_0.01 <- length(AFL_vs_SED_out_genes_0.01)

#Overlap sizes
n12 <- length(intersect(EX_vs_SED_out_genes_0.01, EX_vs_AFL_out_genes_0.01))
n13 <- length(intersect(EX_vs_SED_out_genes_0.01, AFL_vs_SED_out_genes_0.01))
n23 <-  length(intersect(EX_vs_AFL_out_genes_0.01, AFL_vs_SED_out_genes_0.01))
n123 <- length(overlapping_study2_genes)

venn.plot <- draw.triple.venn(
  area1 = length_EX_vs_SED_out_0.01,
  area2 = length_EX_vs_AFL_out_0.01,
  area3 = length_AFL_vs_SED_out_0.01,
  n12 = n12,
  n13 = n13,
  n23 = n23,
  n123 = n123,
  category = c("EXvsSED", "EXvsAFL", "AFLvsSED"),
  fill = c("skyblue", "pink", "salmon"),
  lty = "blank",
  margin = 0.1
)

# Display the Venn diagram
#jpeg("Venndiagram_Study2_genes.jpeg", width = 2500, height = 1800, res = 300)
grid.draw(venn.plot)

# Add a title
grid.text(
  "Venn Diagram for Study2 genes filtered at Pval 0.01", # Title text
  x = 0.5,                                    
  y = 0.95,                                    
  gp = gpar(fontsize = 12, fontface = "bold")  
)
#dev.off()

## code for extracting the venndiagram data
unique_EX_vs_SED_out_genes_0.01 <- setdiff(EX_vs_SED_out_genes_0.01, union(EX_vs_AFL_out_genes_0.01, AFL_vs_SED_out_genes_0.01))
unique_EX_vs_AFL_out_genes_0.01 <- setdiff(EX_vs_AFL_out_genes_0.01, union(EX_vs_SED_out_genes_0.01, AFL_vs_SED_out_genes_0.01))
unique_AFL_vs_SED_out_genes_0.01 <- setdiff(AFL_vs_SED_out_genes_0.01, union(EX_vs_SED_out_genes_0.01, EX_vs_AFL_out_genes_0.01))

# Combine all unique and overlapping genes
all_study2_genes <- unique(c(EX_vs_SED_out_genes_0.01, EX_vs_AFL_out_genes_0.01, AFL_vs_SED_out_genes_0.01))

# Create a data frame for venn data
venn_data_study2 <- data.frame(
  gene_name = all_study2_genes,
  EXvsSED = ifelse(all_study2_genes %in% EX_vs_SED_out_genes_0.01, 1, 0),
  EXvsAFL = ifelse(all_study2_genes %in% EX_vs_AFL_out_genes_0.01, 1, 0),
  AFLvsSED = ifelse(all_study2_genes %in% AFL_vs_SED_out_genes_0.01, 1, 0)
)

#write.xlsx(venn_data_study2, "Venndata_Study2genes.xlsx")

##------------------------------------------------Gene Set Enrichement Analysis - EX vs SED---------------------------------
Hallmark <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt") # from MSigDB
Wiki <- gmtPathways("c2.cp.wikipathways.v2024.1.Hs.symbols.gmt") #from MSigDB
KEGG <- gmtPathways("c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt") #from MSigDB

# prepare ranked gene list for EX vs SED
ranked_genes_EXvsSED <- EX_vs_SED_out$logFC
names(ranked_genes_EXvsSED) <- EX_vs_SED_out$gene_name

# Sort genes in decreasing order of logFC
ranked_genes_EXvsSED <- sort(ranked_genes_EXvsSED, decreasing = TRUE)

# View the ranked genes
head(ranked_genes_EXvsSED)

# Run fgsea with Hallmark pathways 
FGSEA_hl_EXvsSED <- fgsea(pathways = Hallmark , # List of gene sets to check
                      stats = ranked_genes_EXvsSED,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_hl_EXvsSED$leadingEdge <- sapply(FGSEA_hl_EXvsSED$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string
#write.xlsx(FGSEA_hl_EXvsSED, "/Users/JavvadP/Downloads/Wagner/out_fgsea_study2.xlsx")

#select top pathways and create a table
topPathwaysUp_hl_EXvsSED <- FGSEA_hl_EXvsSED[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 hallmark pathways upregulated at pval 0.01
topPathwaysDown_hl_EXvsSED <- FGSEA_hl_EXvsSED[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 hallmark pathways downregulated at pval 0.01
topPathways_hl_EXvsSED <- c(topPathwaysUp_hl_EXvsSED, rev(topPathwaysDown_hl_EXvsSED))
plotGseaTable(Hallmark[topPathways_hl_EXvsSED], stats = ranked_genes_EXvsSED, fgseaRes = FGSEA_hl_EXvsSED, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_hl_EXvsSED <- collapsePathways(FGSEA_hl_EXvsSED[order(pval)][padj < 0.01], Hallmark, ranked_genes_EXvsSED, gseaParam = 1)

# plot the most significantly enriched pathway
plotEnrichment(Hallmark[[head(FGSEA_hl_EXvsSED[order(pval), ], 1)$pathway]],
               ranked_genes_EXvsSED) + 
  labs(title = head(FGSEA_hl_EXvsSED[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_hl_EXvsSED = Hallmark[topPathways_hl_EXvsSED]

# generate data for input into pathway images
df_hl_EXvsSED = gseaCurve(ranked_genes_EXvsSED, setlist_hl_EXvsSED, FGSEA_hl_EXvsSED)

#jpeg("HallmarkPathway_EXvsSED.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_hl_EXvsSED, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "HallmarkPathway_EXvsSED(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()


## Wiki pathways
# Run fgsea with Wiki pathways 
FGSEA_wk_EXvsSED <- fgsea(pathways = Wiki , # List of gene sets to check
                      stats = ranked_genes_EXvsSED,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_wk_EXvsSED$leadingEdge <- sapply(FGSEA_wk_EXvsSED$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string


#select top pathways and create a table
topPathwaysUp_wk_EXvsSED <- FGSEA_wk_EXvsSED[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 wiki pathways upregulated at pval 0.01
topPathwaysDown_wk_EXvsSED <- FGSEA_wk_EXvsSED[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 wiki pathways downregulated at pval 0.01
topPathways_wk_EXvsSED <- c(topPathwaysUp_wk_EXvsSED, rev(topPathwaysDown_wk_EXvsSED))
plotGseaTable(Wiki[topPathways_wk_EXvsSED], stats = ranked_genes_EXvsSED, fgseaRes = FGSEA_wk_EXvsSED, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_wk_EXvsSED <- collapsePathways(FGSEA_wk_EXvsSED[order(pval)][padj < 0.01], Wiki, ranked_genes_EXvsSED, gseaParam = 1)
mainPathways_wk_EXvsSED <- FGSEA_wk_EXvsSED[pathway %in% collapsedPathways_wk_EXvsSED$mainPathways][order(-NES), pathway]
plotGseaTable(Wiki[mainPathways_wk_EXvsSED], ranked_genes_EXvsSED, FGSEA_wk_EXvsSED, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(Wiki[[head(FGSEA_wk_EXvsSED[order(pval), ], 1)$pathway]],
               ranked_genes_EXvsSED) + 
  labs(title = head(FGSEA_wk_EXvsSED[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_wk_EXvsSED = Wiki[topPathways_wk_EXvsSED]

# generate data for input into pathway images
df_wk_EXvsSED = gseaCurve(ranked_genes_EXvsSED, setlist_wk_EXvsSED, FGSEA_wk_EXvsSED)

#jpeg("WikiPathway_EXvsSED.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_wk_EXvsSED, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "WikiPathway_EXvsSED(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

## KEGG pathway
# Run fgsea with KEGG pathways 
FGSEA_kg_EXvsSED <- fgsea(pathways = KEGG , # List of gene sets to check
                      stats = ranked_genes_EXvsSED,
                      eps = 0,
                      scoreType = 'std', 
                      minSize = 15,
                      maxSize = 250)

#FGSEA_kg_EXvsSED$leadingEdge <- sapply(FGSEA_kg_EXvsSED$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string


#select top pathways and create a table
topPathwaysUp_kg_EXvsSED <- FGSEA_kg_EXvsSED[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 kegg pathways upregulated at pval 0.01
topPathwaysDown_kg_EXvsSED <- FGSEA_kg_EXvsSED[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 kegg pathways downregulated at pval 0.01
topPathways_kg_EXvsSED <- c(topPathwaysUp_kg_EXvsSED, rev(topPathwaysDown_kg_EXvsSED))
plotGseaTable(KEGG[topPathways_kg_EXvsSED], stats = ranked_genes_EXvsSED, fgseaRes = FGSEA_kg_EXvsSED, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_kg_EXvsSED <- collapsePathways(FGSEA_kg_EXvsSED[order(pval)][padj < 0.01], KEGG, ranked_genes_EXvsSED, gseaParam = 1)
mainPathways_kg_EXvsSED <- FGSEA_kg_EXvsSED[pathway %in% collapsedPathways_kg_EXvsSED$mainPathways][order(-NES), pathway]
plotGseaTable(KEGG[mainPathways_kg_EXvsSED], ranked_genes_EXvsSED, FGSEA_kg_EXvsSED, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(KEGG[[head(FGSEA_kg_EXvsSED[order(pval), ], 1)$pathway]],
               ranked_genes_EXvsSED) + 
  labs(title = head(FGSEA_kg_EXvsSED[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_kg_EXvsSED = KEGG[topPathways_kg_EXvsSED]

# generate data for input into pathway images
df_kg_EXvsSED = gseaCurve(ranked_genes_EXvsSED, setlist_kg_EXvsSED, FGSEA_kg_EXvsSED)

#jpeg("KEGGPathway_EXvsSED.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_kg_EXvsSED, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "KEGGPathway_EXvsSED(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()


#FGSEA_hl_EXvsSED <- FGSEA_hl_EXvsSED %>% mutate(Database = "Hallmark")
#FGSEA_wk_EXvsSED <- FGSEA_wk_EXvsSED %>% mutate(Database = "Wiki")
#FGSEA_kg_EXvsSED <- FGSEA_kg_EXvsSED %>% mutate(Database = "KEGG")

#combine_EXvsSED <- bind_rows(FGSEA_hl_EXvsSED, FGSEA_wk_EXvsSED, FGSEA_kg_EXvsSED)
#combine_EXvsSED <- combine_EXvsSED %>% select(Database, everything())

#write.xlsx(combine_EXvsSED, "/Users/JavvadP/Downloads/Wagner/out_fgsea_study2.xlsx")

#wb <- loadWorkbook("out_fgsea_study1.xlsx")
#addWorksheet(wb, "TNFvsTCZ_EXvsSED")
#writeData(wb, "TNFvsTCZ_EXvsSED", combine_EXvsSED)
#saveWorkbook(wb, "out_fgsea_Study1.xlsx", overwrite = T)

##--------------------------------------Barplots- FGSEA - EXvsSED------------------------------------------
top_hl_EXvsSED <- FGSEA_hl_EXvsSED[pathway %in% topPathways_hl_EXvsSED][order(pval)] #pull the top pathways from hallmark db

#Barplot for hallmark pathway - EXvsSED

#jpeg("FGSEA_Hallmark_EXvsSED.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_hl_EXvsSED) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Hallmark Pathway - EXvsSED")
#dev.off()


## Wiki pathways
top_wk_EXvsSED <- FGSEA_wk_EXvsSED[pathway %in% topPathways_wk_EXvsSED][order(pval)] #pull the top pathways from wiki db

#Barplot for wiki pathway - EXvsSED

#jpeg("FGSEA_Wiki_EXvsSED.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_wk_EXvsSED) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Wiki Pathway - EXvsSED")
#dev.off()

## KEGG pathways
top_kg_EXvsSED <- FGSEA_kg_EXvsSED[pathway %in% topPathways_kg_EXvsSED][order(pval)] #pull the top pathways from KEGG db

#Barplot for KEGG pathway - EXvsSED

#jpeg("FGSEA_KEGG_EXvsSED.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_kg_EXvsSED) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="KEGG Pathway - EXvsSED")
#dev.off()

##------------------------------------------------Gene Set Enrichement Analysis -EX vs AFL---------------------------------

# prepare ranked gene list for EX vs AFL
ranked_genes_EXvsAFL <- EX_vs_AFL_out$logFC
names(ranked_genes_EXvsAFL) <- EX_vs_AFL_out$gene_name

# Sort genes in decreasing order of logFC
ranked_genes_EXvsAFL <- sort(ranked_genes_EXvsAFL, decreasing = TRUE)

# View the ranked genes
head(ranked_genes_EXvsAFL)

# Run fgsea with Hallmark pathways 
FGSEA_hl_EXvsAFL <- fgsea(pathways = Hallmark , # List of gene sets to check
                          stats = ranked_genes_EXvsAFL,
                          eps = 0,
                          scoreType = 'std', 
                          minSize = 15,
                          maxSize = 250)

#FGSEA_hl_EXvsAFL$leadingEdge <- sapply(FGSEA_hl_EXvsAFL$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string

#select top pathways and create a table
topPathwaysUp_hl_EXvsAFL <- FGSEA_hl_EXvsAFL[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 hallmark pathways upregulated at pval 0.01
topPathwaysDown_hl_EXvsAFL <- FGSEA_hl_EXvsAFL[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 hallmark pathways downregulated at pval 0.01
topPathways_hl_EXvsAFL <- c(topPathwaysUp_hl_EXvsAFL, rev(topPathwaysDown_hl_EXvsAFL))
plotGseaTable(Hallmark[topPathways_hl_EXvsAFL], stats = ranked_genes_EXvsAFL, fgseaRes = FGSEA_hl_EXvsAFL, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_hl_EXvsAFL <- collapsePathways(FGSEA_hl_EXvsAFL[order(pval)][padj < 0.01], Hallmark, ranked_genes_EXvsAFL, gseaParam = 1)

# plot the most significantly enriched pathway
plotEnrichment(Hallmark[[head(FGSEA_hl_EXvsAFL[order(pval), ], 1)$pathway]],
               ranked_genes_EXvsAFL) + 
  labs(title = head(FGSEA_hl_EXvsAFL[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_hl_EXvsAFL = Hallmark[topPathways_hl_EXvsAFL]

# generate data for input into pathway images
df_hl_EXvsAFL = gseaCurve(ranked_genes_EXvsAFL, setlist_hl_EXvsAFL, FGSEA_hl_EXvsAFL)

#jpeg("HallmarkPathway_EXvsAFL.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_hl_EXvsAFL, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "HallmarkPathway_EXvsAFL(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()


## Wiki pathways
# Run fgsea with Wiki pathways 
FGSEA_wk_EXvsAFL <- fgsea(pathways = Wiki , # List of gene sets to check
                          stats = ranked_genes_EXvsAFL,
                          eps = 0,
                          scoreType = 'std', 
                          minSize = 15,
                          maxSize = 250)

#FGSEA_wk_EXvsAFL$leadingEdge <- sapply(FGSEA_wk_EXvsAFL$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string

#select top pathways and create a table
topPathwaysUp_wk_EXvsAFL <- FGSEA_wk_EXvsAFL[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 wiki pathways upregulated at pval 0.01
topPathwaysDown_wk_EXvsAFL <- FGSEA_wk_EXvsAFL[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 wiki pathways downregulated at pval 0.01
topPathways_wk_EXvsAFL <- c(topPathwaysUp_wk_EXvsAFL, rev(topPathwaysDown_wk_EXvsAFL))
plotGseaTable(Wiki[topPathways_wk_EXvsAFL], stats = ranked_genes_EXvsAFL, fgseaRes = FGSEA_wk_EXvsAFL, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_wk_EXvsAFL <- collapsePathways(FGSEA_wk_EXvsAFL[order(pval)][padj < 0.01], Wiki, ranked_genes_EXvsAFL, gseaParam = 1)
mainPathways_wk_EXvsAFL <- FGSEA_wk_EXvsAFL[pathway %in% collapsedPathways_wk_EXvsAFL$mainPathways][order(-NES), pathway]
plotGseaTable(Wiki[mainPathways_wk_EXvsAFL], ranked_genes_EXvsAFL, FGSEA_wk_EXvsAFL, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(Wiki[[head(FGSEA_wk_EXvsAFL[order(pval), ], 1)$pathway]],
               ranked_genes_EXvsAFL) + 
  labs(title = head(FGSEA_wk_EXvsAFL[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_wk_EXvsAFL = Wiki[topPathways_wk_EXvsAFL]

# generate data for input into pathway images
df_wk_EXvsAFL = gseaCurve(ranked_genes_EXvsAFL, setlist_wk_EXvsAFL, FGSEA_wk_EXvsAFL)

#jpeg("WikiPathway_EXvsAFL.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_wk_EXvsAFL, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "WikiPathway_EXvsAFL(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

## KEGG pathway
# Run fgsea with KEGG pathways 
FGSEA_kg_EXvsAFL <- fgsea(pathways = KEGG , # List of gene sets to check
                          stats = ranked_genes_EXvsAFL,
                          eps = 0,
                          scoreType = 'std', 
                          minSize = 15,
                          maxSize = 250)

#FGSEA_kg_EXvsAFL$leadingEdge <- sapply(FGSEA_kg_EXvsAFL$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string

#select top pathways and create a table
topPathwaysUp_kg_EXvsAFL <- FGSEA_kg_EXvsAFL[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 kegg pathways upregulated at pval 0.01
topPathwaysDown_kg_EXvsAFL <- FGSEA_kg_EXvsAFL[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 kegg pathways downregulated at pval 0.01
topPathways_kg_EXvsAFL <- c(topPathwaysUp_kg_EXvsAFL, rev(topPathwaysDown_kg_EXvsAFL))
plotGseaTable(KEGG[topPathways_kg_EXvsAFL], stats = ranked_genes_EXvsAFL, fgseaRes = FGSEA_kg_EXvsAFL, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_kg_EXvsAFL <- collapsePathways(FGSEA_kg_EXvsAFL[order(pval)][padj < 0.01], KEGG, ranked_genes_EXvsAFL, gseaParam = 1)
mainPathways_kg_EXvsAFL <- FGSEA_kg_EXvsAFL[pathway %in% collapsedPathways_kg_EXvsAFL$mainPathways][order(-NES), pathway]
plotGseaTable(KEGG[mainPathways_kg_EXvsAFL], ranked_genes_EXvsAFL, FGSEA_kg_EXvsAFL, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(KEGG[[head(FGSEA_kg_EXvsAFL[order(pval), ], 1)$pathway]],
               ranked_genes_EXvsAFL) + 
  labs(title = head(FGSEA_kg_EXvsAFL[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_kg_EXvsAFL = KEGG[topPathways_kg_EXvsAFL]

# generate data for input into pathway images
df_kg_EXvsAFL = gseaCurve(ranked_genes_EXvsAFL, setlist_kg_EXvsAFL, FGSEA_kg_EXvsAFL)

#jpeg("KEGGPathway_EXvsAFL.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_kg_EXvsAFL, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "KEGGPathway_EXvsAFL(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

#FGSEA_hl_EXvsAFL <- FGSEA_hl_EXvsAFL %>% mutate(Database = "Hallmark")
#FGSEA_wk_EXvsAFL <- FGSEA_wk_EXvsAFL %>% mutate(Database = "Wiki")
#FGSEA_kg_EXvsAFL <- FGSEA_kg_EXvsAFL %>% mutate(Database = "KEGG")

#combine_EXvsAFL <- bind_rows(FGSEA_hl_EXvsAFL, FGSEA_wk_EXvsAFL, FGSEA_kg_EXvsAFL)
#combine_EXvsAFL <- combine_EXvsAFL %>% select(Database, everything())

#wb <- loadWorkbook("/Users/JavvadP/Downloads/Wagner/out_fgsea_study2.xlsx")
#addWorksheet(wb, "EXvsAFL")
#writeData(wb, "EXvsAFL", combine_EXvsAFL)
#saveWorkbook(wb, "/Users/JavvadP/Downloads/Wagner/out_fgsea_study2.xlsx", overwrite = T)

##--------------------------------------Barplots- FGSEA -EXvsAFL------------------------------------------
top_hl_EXvsAFL <- FGSEA_hl_EXvsAFL[pathway %in% topPathways_hl_EXvsAFL][order(pval)] #pull the top pathways from hallmark db

#Barplot for hallmark pathway -EXvsAFL

#jpeg("FGSEA_Hallmark_EXvsAFL.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_hl_EXvsAFL) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Hallmark Pathway -EXvsAFL")
#dev.off()


## Wiki pathways
top_wk_EXvsAFL <- FGSEA_wk_EXvsAFL[pathway %in% topPathways_wk_EXvsAFL][order(pval)] #pull the top pathways from wiki db

#Barplot for wiki pathway -EXvsAFL

#jpeg("FGSEA_Wiki_EXvsAFL.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_wk_EXvsAFL) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Wiki Pathway -EXvsAFL")
#dev.off()

## KEGG pathways
top_kg_EXvsAFL <- FGSEA_kg_EXvsAFL[pathway %in% topPathways_kg_EXvsAFL][order(pval)] #pull the top pathways from KEGG db

#Barplot for KEGG pathway -EXvsAFL

#jpeg("FGSEA_KEGG_EXvsAFL.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_kg_EXvsAFL) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="KEGG Pathway -EXvsAFL")
#dev.off()

##------------------------------------------------Gene Set Enrichement Analysis -AFL vs SED---------------------------------
# prepare ranked gene list for AFL vs SED
ranked_genes_AFLvsSED <- AFL_vs_SED_out$logFC
names(ranked_genes_AFLvsSED) <- AFL_vs_SED_out$gene_name

# Sort genes in decreasing order of logFC
ranked_genes_AFLvsSED <- sort(ranked_genes_AFLvsSED, decreasing = TRUE)

# View the ranked genes
head(ranked_genes_AFLvsSED)

# Run fgsea with Hallmark pathways 
FGSEA_hl_AFLvsSED <- fgsea(pathways = Hallmark , # List of gene sets to check
                          stats = ranked_genes_AFLvsSED,
                          eps = 0,
                          scoreType = 'std', 
                          minSize = 15,
                          maxSize = 250)

#FGSEA_hl_AFLvsSED$leadingEdge <- sapply(FGSEA_hl_AFLvsSED$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string

#select top pathways and create a table
topPathwaysUp_hl_AFLvsSED <- FGSEA_hl_AFLvsSED[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 hallmark pathways upregulated at pval 0.01
topPathwaysDown_hl_AFLvsSED <- FGSEA_hl_AFLvsSED[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 hallmark pathways downregulated at pval 0.01
topPathways_hl_AFLvsSED <- c(topPathwaysUp_hl_AFLvsSED, rev(topPathwaysDown_hl_AFLvsSED))
plotGseaTable(Hallmark[topPathways_hl_AFLvsSED], stats = ranked_genes_AFLvsSED, fgseaRes = FGSEA_hl_AFLvsSED, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_hl_AFLvsSED <- collapsePathways(FGSEA_hl_AFLvsSED[order(pval)][padj < 0.01], Hallmark, ranked_genes_AFLvsSED, gseaParam = 1)

# plot the most significantly enriched pathway
plotEnrichment(Hallmark[[head(FGSEA_hl_AFLvsSED[order(pval), ], 1)$pathway]],
               ranked_genes_AFLvsSED) + 
  labs(title = head(FGSEA_hl_AFLvsSED[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_hl_AFLvsSED = Hallmark[topPathways_hl_AFLvsSED]

# generate data for input into pathway images
df_hl_AFLvsSED = gseaCurve(ranked_genes_AFLvsSED, setlist_hl_AFLvsSED, FGSEA_hl_AFLvsSED)

#jpeg("HallmarkPathway_AFLvsSED.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_hl_AFLvsSED, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "HallmarkPathway_AFLvsSED(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()


## Wiki pathways
# Run fgsea with Wiki pathways 
FGSEA_wk_AFLvsSED <- fgsea(pathways = Wiki , # List of gene sets to check
                          stats = ranked_genes_AFLvsSED,
                          eps = 0,
                          scoreType = 'std', 
                          minSize = 15,
                          maxSize = 250)

#FGSEA_wk_AFLvsSED$leadingEdge <- sapply(FGSEA_wk_AFLvsSED$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string

#select top pathways and create a table
topPathwaysUp_wk_AFLvsSED <- FGSEA_wk_AFLvsSED[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 wiki pathways upregulated at pval 0.01
topPathwaysDown_wk_AFLvsSED <- FGSEA_wk_AFLvsSED[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 wiki pathways downregulated at pval 0.01
topPathways_wk_AFLvsSED <- c(topPathwaysUp_wk_AFLvsSED, rev(topPathwaysDown_wk_AFLvsSED))
plotGseaTable(Wiki[topPathways_wk_AFLvsSED], stats = ranked_genes_AFLvsSED, fgseaRes = FGSEA_wk_AFLvsSED, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_wk_AFLvsSED <- collapsePathways(FGSEA_wk_AFLvsSED[order(pval)][padj < 0.01], Wiki, ranked_genes_AFLvsSED, gseaParam = 1)
mainPathways_wk_AFLvsSED <- FGSEA_wk_AFLvsSED[pathway %in% collapsedPathways_wk_AFLvsSED$mainPathways][order(-NES), pathway]
plotGseaTable(Wiki[mainPathways_wk_AFLvsSED], ranked_genes_AFLvsSED, FGSEA_wk_AFLvsSED, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(Wiki[[head(FGSEA_wk_AFLvsSED[order(pval), ], 1)$pathway]],
               ranked_genes_AFLvsSED) + 
  labs(title = head(FGSEA_wk_AFLvsSED[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_wk_AFLvsSED = Wiki[topPathways_wk_AFLvsSED]

# generate data for input into pathway images
df_wk_AFLvsSED = gseaCurve(ranked_genes_AFLvsSED, setlist_wk_AFLvsSED, FGSEA_wk_AFLvsSED)

#jpeg("WikiPathway_AFLvsSED.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_wk_AFLvsSED, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "WikiPathway_AFLvsSED(P.val<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

## KEGG pathway
# Run fgsea with KEGG pathways 
FGSEA_kg_AFLvsSED <- fgsea(pathways = KEGG , # List of gene sets to check
                          stats = ranked_genes_AFLvsSED,
                          eps = 0,
                          scoreType = 'std', 
                          minSize = 15,
                          maxSize = 250)

#FGSEA_kg_AFLvsSED$leadingEdge <- sapply(FGSEA_kg_AFLvsSED$leadingEdge, function(x) paste(x, collapse=",")) #converting each element of Leading Edge into character string

#select top pathways and create a table
topPathwaysUp_kg_AFLvsSED <- FGSEA_kg_AFLvsSED[ES > 0 & pval < 0.01][head(order(pval), n = 5), pathway] # Extract the top 5 kegg pathways upregulated at pval 0.01
topPathwaysDown_kg_AFLvsSED <- FGSEA_kg_AFLvsSED[ES < 0 & pval < 0.01][head(order(pval), n =5), pathway] # Extract the top 5 kegg pathways downregulated at pval 0.01
topPathways_kg_AFLvsSED <- c(topPathwaysUp_kg_AFLvsSED, rev(topPathwaysDown_kg_AFLvsSED))
plotGseaTable(KEGG[topPathways_kg_AFLvsSED], stats = ranked_genes_AFLvsSED, fgseaRes = FGSEA_kg_AFLvsSED, gseaParam = 0.5)

# remove redundant pathways
collapsedPathways_kg_AFLvsSED <- collapsePathways(FGSEA_kg_AFLvsSED[order(pval)][padj < 0.01], KEGG, ranked_genes_AFLvsSED, gseaParam = 1)
mainPathways_kg_AFLvsSED <- FGSEA_kg_AFLvsSED[pathway %in% collapsedPathways_kg_AFLvsSED$mainPathways][order(-NES), pathway]
plotGseaTable(KEGG[mainPathways_kg_AFLvsSED], ranked_genes_AFLvsSED, FGSEA_kg_AFLvsSED, gseaParam = 0.5)


# plot the most significantly enriched pathway
plotEnrichment(KEGG[[head(FGSEA_kg_AFLvsSED[order(pval), ], 1)$pathway]],
               ranked_genes_AFLvsSED) + 
  labs(title = head(FGSEA_kg_AFLvsSED[order(pval), ], 1)$pathway)

# using gggsea to generate enrichment plots for top pathways (identified from fgsea above)

# call selected pathways file
setlist_kg_AFLvsSED = KEGG[topPathways_kg_AFLvsSED]

# generate data for input into pathway images
df_kg_AFLvsSED = gseaCurve(ranked_genes_AFLvsSED, setlist_kg_AFLvsSED, FGSEA_kg_AFLvsSED)

#jpeg("KEGGPathway_AFLvsSED.jpeg", width = 2500, height = 1800, res = 300)
ggplot2::ggplot() +
  geom_gsea(df_kg_AFLvsSED, linecolor = "darkgreen", zeroline = F, linesize = 1, labelsize = 2) +
  theme_gsea(7) +
  labs(title = "KEGGPathway_AFLvsSED(Pval<=0.01)") +
  theme(plot.title = element_text(face = "bold", size = 13))
#dev.off()

#FGSEA_hl_AFLvsSED <- FGSEA_hl_AFLvsSED %>% mutate(Database = "Hallmark")
#FGSEA_wk_AFLvsSED <- FGSEA_wk_AFLvsSED %>% mutate(Database = "Wiki")
#FGSEA_kg_AFLvsSED <- FGSEA_kg_AFLvsSED %>% mutate(Database = "KEGG")

#combine_AFLvsSED <- bind_rows(FGSEA_hl_AFLvsSED, FGSEA_wk_AFLvsSED, FGSEA_kg_AFLvsSED)
#combine_AFLvsSED <- combine_AFLvsSED %>% select(Database, everything())

#wb <- loadWorkbook("/Users/JavvadP/Downloads/Wagner/out_fgsea_study2.xlsx")
#addWorksheet(wb, "AFLvsSED")
#writeData(wb, "AFLvsSED", combine_AFLvsSED)
#saveWorkbook(wb, "/Users/JavvadP/Downloads/Wagner/out_fgsea_study2.xlsx", overwrite = T)


##--------------------------------------Barplots- FGSEA -AFLvsSED------------------------------------------
top_hl_AFLvsSED <- FGSEA_hl_AFLvsSED[pathway %in% topPathways_hl_AFLvsSED][order(pval)] #pull the top pathways from hallmark db

#Barplot for hallmark pathway -AFLvsSED

#jpeg("FGSEA_Hallmark_AFLvsSED.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_hl_AFLvsSED) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Hallmark Pathway -AFLvsSED")
#dev.off()


## Wiki pathways
top_wk_AFLvsSED <- FGSEA_wk_AFLvsSED[pathway %in% topPathways_wk_AFLvsSED][order(pval)] #pull the top pathways from wiki db

#Barplot for wiki pathway -AFLvsSED

#jpeg("FGSEA_Wiki_AFLvsSED.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_wk_AFLvsSED) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="Wiki Pathway -AFLvsSED")
#dev.off()

## KEGG pathways
top_kg_AFLvsSED <- FGSEA_kg_AFLvsSED[pathway %in% topPathways_kg_AFLvsSED][order(pval)] #pull the top pathways from KEGG db

#Barplot for KEGG pathway -AFLvsSED

#jpeg("FGSEA_KEGG_AFLvsSED.jpeg", width = 2500, height = 1800, res = 300)
ggplot(top_kg_AFLvsSED) +
  geom_col(aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = -log10(pval)
  ), color = "black") +
  coord_flip() +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.line.x = element_line(color="gray30", linetype = "solid")
  ) +
  geom_hline (yintercept = 0, color="black", linetype = "solid")+
  labs(y="Pathway Regulation (NES)", x="Pathway Name", 
       title="KEGG Pathway -AFLvsSED")
#dev.off()


##-----------------------------------------New analysis - 2/10/25 Pathway overlap of study2 comparisons-------------------
# get the EXvsSED unique pathways
Hallmark_EXvsSED <- unique(topPathways_hl_EXvsSED)
Wiki_EXvsSED <- unique(topPathways_wk_EXvsSED)
Kegg_EXvsSED <- unique(topPathways_kg_EXvsSED)

#get the EXvsAFL unique pathways
Hallmark_EXvsAFL <- unique(topPathways_hl_EXvsAFL)
Wiki_EXvsAFL <- unique(topPathways_wk_EXvsAFL)
Kegg_EXvsAFL <- unique(topPathways_kg_EXvsAFL)

#get the AFlvsSED unique pathways
Hallmark_AFLvsSED <- unique(topPathways_hl_AFLvsSED)
wiki_AFLvsSED <- unique(topPathways_wk_AFLvsSED)
Kegg_AFLvsSED <- unique(topPathways_kg_AFLvsSED)

hl_list <- list(
  Hallmark_EXvsSED = Hallmark_EXvsSED,
  Hallmark_EXvsAFL = Hallmark_EXvsAFL,
  Hallmark_AFLvsSED = Hallmark_AFLvsSED
)

# Generate venn Diagram
#jpeg("Venndiagram_HallmarkPathway_Study2.jpeg", width = 2500, height = 1800, res = 300)
ggVennDiagram(hl_list, category.names = c("EXvsSED", "EXvsAFL", "AFLvsSED"),
              label = "count",set_color = c("blue","red","darkgreen"),
              label_color = "black", label_size = 4,) +
  scale_x_continuous(expand = expansion(mult =.4)) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  labs(title   = "HallmarkPathway overlap between study2 comparisons 
                (EXvsSED, EXvsAFL, AFLvsSED)")
#dev.off()

wk_list <- list(
  Wiki_EXvsSED = Wiki_EXvsSED,
  Wiki_EXvsAFL = Wiki_EXvsAFL,
  wiki_AFLvsSED = wiki_AFLvsSED
)


# Generate venn Diagram
#jpeg("Venndiagram_WikiPathway_Study2.jpeg", width = 2500, height = 1800, res = 300)
ggVennDiagram(wk_list, category.names = c("EXvsSED", "EXvsAFL", "AFLvsSED"),
              label = "count",set_color = c("blue","red","darkgreen"),
              label_color = "black", label_size = 4,) +
  scale_x_continuous(expand = expansion(mult =.4)) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  labs(title   = "WikiPathway overlap between study2 comparisons 
                (EXvsSED, EXvsAFL, AFLvsSED)")
#dev.off()

kg_list<- list(
  Kegg_EXvsSED = Kegg_EXvsSED,
  Kegg_EXvsAFL = Kegg_EXvsAFL,
  Kegg_AFLvsSED = Kegg_AFLvsSED
)

# Generate venn Diagram
#jpeg("Venndiagram_KEGGPathway_Study2.jpeg", width = 2500, height = 1800, res = 300)
ggVennDiagram(wk_list, category.names = c("EXvsSED", "EXvsAFL", "AFLvsSED"),
              label = "count",set_color = c("blue","red","darkgreen"),
              label_color = "black", label_size = 4,) +
  scale_x_continuous(expand = expansion(mult =.4)) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  labs(title   = "KEGGPathway overlap between study2 comparisons 
                (EXvsSED, EXvsAFL, AFLvsSED)")
#dev.off()


## Extract the venn data
Hallmark_EXvsSED <- data.frame(Pathway = Hallmark_EXvsSED)
Hallmark_EXvsSED <- Hallmark_EXvsSED %>% select(Pathway) %>% mutate(Hallmark_EXvsSED = 1)
Hallmark_EXvsAFL <- data.frame(Pathway = Hallmark_EXvsAFL)
Hallmark_EXvsAFL <- Hallmark_EXvsAFL %>% select(Pathway) %>% mutate(Hallmark_EXvsAFL = 1)
Hallmark_AFLvsSED <- data.frame(Pathway = Hallmark_AFLvsSED)
Hallmark_AFLvsSED <- Hallmark_AFLvsSED %>% select(Pathway) %>% mutate(Hallmark_AFLvsSED = 1)

hl_combine <- list(Hallmark_EXvsSED, Hallmark_EXvsAFL, Hallmark_AFLvsSED) %>%
  reduce(full_join, by = "Pathway") %>%
  replace(is.na(.), 0) # Replace NA values with 0

#write.xlsx(hl_combine, "venndata_study2_pathways.xlsx")


Wiki_EXvsSED <- data.frame(Pathway = Wiki_EXvsSED)
Wiki_EXvsSED <- Wiki_EXvsSED %>% select(Pathway) %>% mutate(Wiki_EXvsSED = 1)
Wiki_EXvsAFL <- data.frame(Pathway = Wiki_EXvsAFL)
Wiki_EXvsAFL <- Wiki_EXvsAFL %>% select(Pathway) %>% mutate(Wiki_EXvsAFL = 1)
wiki_AFLvsSED <- data.frame(Pathway = wiki_AFLvsSED)
wiki_AFLvsSED <- wiki_AFLvsSED %>% select(Pathway) %>% mutate(wiki_AFLvsSED = 1)

wk_combine <- list(Wiki_EXvsSED, Wiki_EXvsAFL, wiki_AFLvsSED) %>%
  reduce(full_join, by = "Pathway") %>%
  replace(is.na(.), 0) # Replace NA values with 0


Kegg_EXvsSED <- data.frame(Pathway = Kegg_EXvsSED)
Kegg_EXvsSED <- Kegg_EXvsSED %>% select(Pathway) %>% mutate(Kegg_EXvsSED = 1)
Kegg_EXvsAFL <- data.frame(Pathway = Kegg_EXvsAFL)
Kegg_EXvsAFL <- Kegg_EXvsAFL %>% select(Pathway) %>% mutate(Kegg_EXvsAFL = 1)
Kegg_AFLvsSED <- data.frame(Pathway = Kegg_AFLvsSED)
Kegg_AFLvsSED <- Kegg_AFLvsSED %>% select(Pathway) %>% mutate(Kegg_AFLvsSED = 1)

kg_combine <- list(Kegg_EXvsSED, Kegg_EXvsAFL, Kegg_AFLvsSED) %>%
  reduce(full_join, by = "Pathway") %>%
  replace(is.na(.), 0) # Replace NA values with 0

wb <- loadWorkbook("venndata_study2_pathways.xlsx")
addWorksheet(wb, "keggpathways")
writeData(wb, "keggpathways", kg_combine)
saveWorkbook(wb, "venndata_study2_pathways.xlsx", overwrite = T)

##------------------------------------------------Overlap of Study1 and Study2 pathways--------------------------
# get the TNF unique pathways 
TNF_Hallmark <- unique(topPathways_hl_TNF)
TNF_Wiki <- unique(topPathways_wk_TNF)
TNF_KEGG <- unique(topPathways_kg_TNF)

# get the TCZ unique pathways 
TCZ_Hallmark <- unique(topPathways_hl_TCZ)
TCZ_Wiki <- unique(topPathways_wk_TCZ)
TCZ_KEGG <- unique(topPathways_kg_TCZ)

# get the EXvsSED unique pathways 
EXvsSED_Hallmark <- unique(topPathways_hl_EXvsSED)
EXvsSED_Wiki <- unique(topPathways_wk_EXvsSED)
EXvsSED_KEGG <- unique(topPathways_kg_EXvsSED)

# get the EXvsAFL unique pathways
EXvsAFL_Hallmark <- unique(topPathways_hl_EXvsAFL)
EXvsAFL_Wiki <- unique(topPathways_wk_EXvsAFL)
EXvsAFL_KEGG <- unique(topPathways_kg_EXvsAFL)

# get the AFLvsSED unique pathways
AFLvsSED_Hallmark <- unique(topPathways_hl_AFLvsSED)
AFLvsSED_Wiki <- unique(topPathways_wk_AFLvsSED)
AFLvsSED_KEGG <- unique(topPathways_kg_AFLvsSED)

## Get the venn data for overlap between study1 and study2 Hallmark pathways
TNF_Hallmark <- data.frame(Pathway = TNF_Hallmark)
TNF_Hallmark <- TNF_Hallmark %>% select(Pathway) %>% mutate(Hallmark_TNF = 1)
TCZ_Hallmark <- data.frame(Pathway = TCZ_Hallmark)
TCZ_Hallmark <- TCZ_Hallmark %>% select(Pathway) %>% mutate(Hallmark_TCZ = 1)
EXvsSED_Hallmark <- data.frame(Pathway = EXvsSED_Hallmark)
EXvsSED_Hallmark <- EXvsSED_Hallmark %>% select(Pathway) %>% mutate(Hallmark_EXvsSED = 1)
EXvsAFL_Hallmark <- data.frame(Pathway = EXvsAFL_Hallmark)
EXvsAFL_Hallmark <- EXvsAFL_Hallmark %>% select(Pathway) %>% mutate(Hallmark_EXvsAFL = 1)
AFLvsSED_Hallmark <- data.frame(Pathway = AFLvsSED_Hallmark)
AFLvsSED_Hallmark <- AFLvsSED_Hallmark %>% select(Pathway) %>% mutate(Hallmark_AFLvsSED = 1)

combine_hl <- list(TNF_Hallmark, TCZ_Hallmark, EXvsSED_Hallmark, EXvsAFL_Hallmark, AFLvsSED_Hallmark) %>%
  reduce(full_join, by = "Pathway") %>%
  replace(is.na(.), 0) # Replace NA values with 0

#jpeg("VennDiagram_HallmarkPathway.jpeg", width = 2500, height = 1800, res = 300)
vennDiagram(combine_hl[, 2:6], names = c("Hallmark_TNF", "Hallmark_TCZ", "Hallmark_EXvsSED", "Hallmark_EXvsAFL", "Hallmark_AFLvsSED"), circle.col=c("red", "green", "blue", "violet", "orange"), cex = 0.7 )
#dev.off()

#write.xlsx(combine_hl, "Overlapping_FGSEA_pathways.xlsx")

## Get the venn data for overlap between study1 and study2 Wiki pathways
## Wiki pathways 
TNF_Wiki <- data.frame(Pathway = TNF_Wiki)
TNF_Wiki <- TNF_Wiki %>% select(Pathway) %>% mutate(Wiki_TNF = 1)
TCZ_Wiki <- data.frame(Pathway = TCZ_Wiki)
TCZ_Wiki <- TCZ_Wiki %>% select(Pathway) %>% mutate(Wiki_TCZ = 1)
EXvsSED_Wiki <- data.frame(Pathway = EXvsSED_Wiki)
EXvsSED_Wiki <- EXvsSED_Wiki %>% select(Pathway) %>% mutate(Wiki_EXvsSED = 1)
EXvsAFL_Wiki <- data.frame(Pathway = EXvsAFL_Wiki)
EXvsAFL_Wiki <- EXvsAFL_Wiki %>% select(Pathway) %>% mutate(Wiki_EXvsAFL = 1)
AFLvsSED_Wiki <- data.frame(Pathway = AFLvsSED_Wiki)
AFLvsSED_Wiki <- AFLvsSED_Wiki %>% select(Pathway) %>% mutate(Wiki_AFLvsSED = 1)


combine_wk <- list(TNF_Wiki, TCZ_Wiki, EXvsSED_Wiki, EXvsAFL_Wiki, AFLvsSED_Wiki) %>%
  reduce(full_join, by = "Pathway") %>%
  replace(is.na(.), 0) # Replace NA values with 0

#jpeg("VennDiagram_WikiPathway.jpeg", width = 2500, height = 1800, res = 300)
vennDiagram(combine_wk[, 2:6], names = c("Wiki_TNF", "Wiki_TCZ", "Wiki_EXvsSED", "Wiki_EXvsAFL", "Wiki_AFLvsSED"), circle.col=c("red", "green", "blue", "violet", "orange"), cex = 0.7 )
#dev.off()

## Get the venn data for overlap between study1 and study2 KEGG pathways
TNF_KEGG <- data.frame(Pathway = TNF_KEGG)
TNF_KEGG <- TNF_KEGG %>% select(Pathway) %>% mutate(KEGG_TNF = 1)
TCZ_KEGG <- data.frame(Pathway = TCZ_KEGG)
TCZ_KEGG <- TCZ_KEGG %>% select(Pathway) %>% mutate(KEGG_TCZ = 1)
EXvsSED_KEGG <- data.frame(Pathway = EXvsSED_KEGG)
EXvsSED_KEGG <- EXvsSED_KEGG %>% select(Pathway) %>% mutate(KEGG_EXvsSED = 1)
EXvsAFL_KEGG <- data.frame(Pathway = EXvsAFL_KEGG)
EXvsAFL_KEGG <- EXvsAFL_KEGG %>% select(Pathway) %>% mutate(KEGG_EXvsAFL = 1)
AFLvsSED_KEGG <- data.frame(Pathway = AFLvsSED_KEGG)
AFLvsSED_KEGG <- AFLvsSED_KEGG %>% select(Pathway) %>% mutate(KEGG_AFLvsSED = 1)

combine_kg <- list(TNF_KEGG, TCZ_KEGG, EXvsSED_KEGG, EXvsAFL_KEGG, AFLvsSED_KEGG) %>%
  reduce(full_join, by = "Pathway") %>%
  replace(is.na(.), 0) # Replace NA values with 0

#jpeg("VennDiagram_KEGGPathway.jpeg", width = 2500, height = 1800, res = 300)
vennDiagram(combine_kg[, 2:6], names = c("KEGG_TNF", "KEGG_TCZ", "KEGG_EXvsSED", "KEGG_EXvsAFL", "KEGG_AFLvsSED"), circle.col=c("red", "green", "blue", "violet", "orange"), cex = 0.7 )
#dev.off()

#wb_pathwys <- loadWorkbook("Overlapping_FGSEA_pathways.xlsx")
#addWorksheet(wb_pathwys, "KEGGPathways")
#writeData(wb_pathwys, "KEGGPathways", combine_kg)
#saveWorkbook(wb_pathwys, "overlapping_FGSEA_pathways.xlsx", overwrite = T)

##----------------------------------------overlapped genes study1 and study2---------------------------
study1_genes <- data.frame(gene_name = c(TNF_PostvsPre$gene_name, TCZ_PostvsPre$gene_name))
study2_genes <- data.frame(gene_name = c(EX_vs_SED_out$gene_name, EX_vs_AFL_out$gene_name, AFL_vs_SED_out$gene_name))

# Identify overlapping genes
overlapping_genes <- intersect(study1_genes, study2_genes)

# Number of genes in each set and overlap
length_study1 <- length(study1_genes$gene_name)
length_study2 <- length(study2_genes$gene_name)
length_overlap_study1_study2_genes <- length(overlapping_genes$gene_name)

venn.plot <- draw.pairwise.venn(
  area1 = length_study1,
  area2 = length_study2,
  cross.area = length_overlap_study1_study2_genes,
  category = c("study1", "study2"),
  fill = c("skyblue", "pink"),
  lty = "blank"
)

##--------------------------------------Gene level attempt to find out if the two studies are linked------------------
# Combine all comparisons into a list
all_comparisons <- list(TNF = TNF_PostvsPre, TCZ = TCZ_PostvsPre, 
                        EXvsSED = EX_vs_SED_out, EXvsAFL = EX_vs_AFL_out, 
                        AFLvsSED = AFL_vs_SED_out)

# filter genes at pval 0.001 from all comparisons
significant_genes <- unique(unlist(lapply(all_comparisons, function(df) {
  df %>% filter(P.Value < 0.001) %>% pull(gene_name)
})))



# Create a data frame with gene names and their maximum p-value across comparisons
gene_pval_summary <- data.frame(
  gene_name = significant_genes,
  highest_pval = sapply(significant_genes, function(gene) {
    pvals <- unlist(lapply(all_comparisons, function(df) {  # Extract all p-values for the gene from all comparisons
      if (gene %in% df$gene_name) {
        df %>% filter(gene_name == gene) %>% pull(P.Value)
      } else {
        NA  # If the gene is not found, return NA
      }
    }))
    min(pvals, na.rm = TRUE)  # Take the min, ignoring NA values
  })
)



# Sort the genes by their min p-value in ascending order
sorted_genes <- gene_pval_summary %>% arrange(highest_pval)

# Extract the top 100 genes
top_100_genes <- head(sorted_genes, 100)

# Extract logFC values for significant genes across all comparisons
logfc_matrix <- sapply(all_comparisons, function(df) {
  # Match significant genes and extract logFC
  match_order <- match(sorted_genes$gene_name, df$gene_name)
  df$logFC[match_order]
})

# Assign row names as significant genes
rownames(logfc_matrix) <- sorted_genes$gene_name

logfc_matrix[is.na(logfc_matrix)] <- NA

logfc_matrix <- na.omit(logfc_matrix)

logfc_matrix_100 <- head(logfc_matrix, 100)


# Compute the distance matrix
dist_genes <- dist(logfc_matrix, method = "euclidean")
dist_genes_100 <- dist(logfc_matrix_100, method = "euclidean")

# Perform hierarchical clustering using Ward.D2 -clustering is based on min the variation within the clusters
hclust_genes <- hclust(dist_genes, method = "ward.D2")
hclust_genes_100 <- hclust(dist_genes_100, method = "ward.D2")

# Plot the dendrogram
# Cut tree into clusters 
#clusters <- cutree(hclust_genes, k = 10)

# Add cluster labels to the data for analysis
#gene_clusters <- data.frame(Gene = rownames(logfc_matrix), Cluster = clusters)

# Plot the dendrogram with colored branches
#library(dendextend)
#dend <- as.dendrogram(hclust_genes)
#dend <- color_branches(dend, k = 10)
#plot(dend, main = "Hierarchical Clustering with Clusters", xlab = "", sub = "")


#jpeg("HierarchialClustering_allgenes_Pval0.001.jpeg", width = 2500, height = 1800, res = 300)
plot(hclust_genes, main = "Hierarchical Clustering of genes at Pval 0.001", xlab = "", sub = "")
#dev.off()

#jpeg("HierarchialClustering_top100_Pval0.001.jpeg", width = 2500, height = 1800, res = 300)
plot(hclust_genes_100, main = "Hierarchial Clustering of top 100 genes at Pval 0.001", xlab = "", sub = "")
rect.hclust(hclust_genes_100, k = 4, border = "red")
#dev.off()

#library(gplots)

#jpeg("HeirarchialClustering_Pval0.001.jpeg", width = 3500, height = 1800, res = 300)
#heatmap.2(
#  logfc_matrix,
#  Rowv = as.dendrogram(hclust_genes),  
#  Colv = FALSE,                       
#  scale = "none",                      
#  trace = "none",                      # No trace lines
#  col = colorRampPalette(c("blue", "white", "red"))(50),
#  breaks = seq(-5.231105, 7.992879, length.out= 51),  #full range
#  main = "Heatmap of genes at Pval 0.001 Across Comparisons",
#  xlab = "Comparisons",
#  ylab = "Genes",
#  cexCol = 1,
#  margins = c(6,10)
#)
#dev.off()

#jpeg("HeirarchialClustering_top100_Pval0.001.jpeg", width = 3800, height = 1800, res = 300)
#heatmap.2(
#  logfc_matrix_100,
#  Rowv = as.dendrogram(hclust_genes_100),  
#  Colv = FALSE,                       
#  scale = "none",                       
#  trace = "none",                      # No trace lines
#  col = colorRampPalette(c("blue", "white", "red"))(50),
#  breaks = seq(-2.834076, 7.992879, length.out= 51),  #full range
#  main = "Heatmap of top 100 genes at Pval 0.001 Across Comparisons",
#  xlab = "Comparisons",
#  ylab = "Genes",
#  cexCol = 1,
#  margins = c(6,10)
#)
#dev.off()


# Define breaks to ensure 0 maps to white
breaks <- c(seq(-5, -0.01, length.out = 25),
            0,
            seq(0.01, 5, length.out = 25))



library(pheatmap)
#jpeg("S1_S2_HC_Pval0.001.jpeg", width = 3500, height = 2500, res = 300)
# Generate the heatmap
pheatmap(
  logfc_matrix,                  # Input matrix
  cluster_rows = hclust_genes,          # Use the dendrogram for rows
  cluster_cols = F,                
  color = colorRampPalette(c("blue", "white", "red"))(50),  # Color gradient
  border_color = NA,
  breaks = breaks,
  main = "Heatmap of genes at Pval 0.001 Across Comparisons", 
  fontsize_col = 10,                   
  fontsize_row = 6,                    
  labels_col = colnames(logfc_matrix),  
  labels_row = rownames(logfc_matrix), 
  legend = TRUE,                
  angle_col = 45,
  treeheight_row = 200
)
     #dev.off()

# Define breaks to ensure 0 maps to white
breaks_100 <- c(seq(-4, -0.01, length.out = 25), 0, seq(0.01, 4, length.out = 25))

#jpeg("S1_S2_HC_100_Pval0.001.jpeg", width = 3500, height = 2500, res = 300)
# Generate the heatmap
pheatmap(
  logfc_matrix_100,                  # Input matrix
  cluster_rows = hclust_genes_100,          # Use the dendrogram for rows
  cluster_cols = FALSE,                
  color = colorRampPalette(c("blue", "white", "red"))(50),  # Color gradient
  border_color = NA,
  breaks = breaks_100,
  main = "Heatmap of top 100 genes at Pval 0.001 Across study1 and study2 
          Comparisons", 
  fontsize_col = 10,                   
  fontsize_row = 6,                    
  labels_col = colnames(logfc_matrix_100),  
  labels_row = rownames(logfc_matrix_100), 
  legend = TRUE,                
  angle_col = 45,
  treeheight_row = 200
)
#dev.off()

##PCA on the 1700 genes
#study1_TNF <- TNF_logCPM_file[TNF_logCPM_file$Gene %in% rownames(logfc_matrix_100), ]
#rownames(study1_TNF) <- NULL
#study1_TCZ <- TCZ_logCPM_file[TCZ_logCPM_file$Gene %in% rownames(logfc_matrix_100), ]
#rownames(study1_TCZ) <- NULL
#study2 <- study2_logCPM_out_file[study2_logCPM_out_file$Gene %in% rownames(logfc_matrix_100), ]
#rownames(study2) <- NULL

#combined_study <- study1_TNF %>%
#  full_join(study1_TCZ, by = "Gene") %>%
#  full_join(study2, by = "Gene")

#write.xlsx(combined_study, "combined_logcpm_100.xlsx")

#TNF_0.05 <- TNF_PostvsPre %>%
#  filter(P.Value <= 0.05) %>%  # Filter genes with p-value ≤ 0.05
#  arrange(P.Value) %>%           # Sort by p-value 
#  slice_head(n = 100)

#TCZ_0.05 <- TCZ_PostvsPre %>%
#  filter(P.Value <= 0.05) %>%
#  arrange(P.Value) %>%
#  slice_head( n = 100)

#EX_vs_SED_0.05 <- EX_vs_SED_out %>%
#  filter(P.Value <= 0.05) %>%
#  arrange(P.Value) %>%
#  slice_head(n = 100)

#EX_vs_AFL_0.05 <- EX_vs_AFL_out %>%
#  filter(P.Value <= 0.05) %>%
#  arrange(P.Value) %>%
#  slice_head(n = 100)

#AFL_vs_SED_0.05 <- AFL_vs_SED_out %>%
#  filter(P.Value <= 0.05) %>%
#  arrange(P.Value) %>%
#  slice_head(n = 100)

#S1_S2_genes <- data.frame(gene_name = c(TNF_0.05$gene_name, TCZ_0.05$gene_name, 
#                                        EX_vs_SED_0.05$gene_name, EX_vs_AFL_0.05$gene_name, 
#                                        AFL_vs_SED_0.05$gene_name))
#S1_S2_genes <- unique(S1_S2_genes)

#all_comparisons_0.05 <- list(TNF = TNF_0.05, TCZ = TCZ_0.05, EX_vs_SED = EX_vs_SED_0.05, EX_vs_AFL = EX_vs_AFL_0.05, AFL_vs_SED = AFL_vs_SED_0.05)

# Extract logFC values for significant genes across all comparisons
#logfc_matrix_0.05 <- sapply(all_comparisons_0.05, function(df) {
  # Match significant genes and extract logFC
#  match_order <- match(S1_S2_genes$gene_name, df$gene_name)
#  df$logFC[match_order]
#})

# Assign row names as significant genes
#rownames(logfc_matrix_0.05) <- S1_S2_genes$gene_name


#filtered_genes <- list(
#  TNF_PostvsPre %>% filter(P.Value <= 0.2) %>% arrange(P.Value) %>% pull(gene_name),
#  TCZ_PostvsPre %>% filter(P.Value <= 0.2) %>% arrange(P.Value) %>% pull(gene_name),
#  EX_vs_SED_out %>% filter(P.Value <= 0.2) %>% arrange(P.Value) %>% pull(gene_name),
#  EX_vs_AFL_out %>% filter(P.Value <= 0.2) %>% arrange(P.Value) %>% pull(gene_name),
#  AFL_vs_SED_out %>% filter(P.Value <= 0.2) %>% arrange(P.Value) %>% pull(gene_name)
#)

# Find genes that are common to all 5 comparisons
#common_genes <- Reduce(intersect, filtered_genes)

#common_genes <- unique(common_genes)

#common_genes <- data.frame(gene_name = common_genes)


#all_comparisons_0.2 <- list(
# TNF =  TNF_PostvsPre %>% filter(P.Value <= 0.2) %>% arrange(P.Value), 
#  TCZ = TCZ_PostvsPre %>% filter(P.Value <= 0.2) %>% arrange(P.Value), 
#  EXvsSED = EX_vs_SED_out %>% filter(P.Value <= 0.2) %>% arrange(P.Value), 
#  EXvsAFL = EX_vs_AFL_out %>% filter(P.Value <= 0.2) %>% arrange(P.Value), 
#  AFLvsSED = AFL_vs_SED_out %>% filter(P.Value <= 0.2) %>% arrange(P.Value))


#common_genes_out <- lapply(all_comparisons_0.2, function(df) {
#       df %>% filter(gene_name %in% common_genes$gene_name)
#})

#top5_genes <- data.frame(gene_name = common_genes_out$EXvsAFL %>%
#  arrange(P.Value) %>%
#  slice_head(n = 5) %>%
#  pull(gene_name))

#logfc_matrix_0.2 <- sapply(all_comparisons_0.2, function(df) {
  # Match significant genes and extract logFC
#    match_order <- match(top5_genes$gene_name, df$gene_name)
#    df$logFC[match_order]
#})

#rownames(logfc_matrix_0.2) <- top5_genes$gene_name

#jpeg("Heatmap_pval0.01.jpeg", width = 2500, height = 1800, res = 300)
#pheatmap(
#  logfc_matrix_0.2,                  # Input matrix
#  cluster_rows = T,         # Use the dendrogram for rows
#  cluster_cols = F,                
#  color = colorRampPalette(c("blue", "white", "red"))(50),  # Color gradient
#  border_color = NA,
#  breaks = seq(-1,1,length = 51),
#  main = "Heatmap of genes at Pval 0.01 Across Comparisons", 
#  fontsize_col = 10,                   
#  fontsize_row = 6,                    
#  labels_col = colnames(logfc_matrix_0.2),  
#  labels_row = rownames(logfc_matrix_0.2), 
#  legend = TRUE,                
#  angle_col = 45
#)
#dev.off()

##-----------------------------------------------------------study2- Analysis2-----------------------------------------
study2_comparisons <- list(EXvsSED = EX_vs_SED_out, EXvsAFL = EX_vs_AFL_out, AFLvsSED = AFL_vs_SED_out)

# filter genes at pval 0.01 from study2 comparisons
study2_significant_genes <- unique(unlist(lapply(study2_comparisons, function(df) {
  df %>% filter(P.Value < 0.01) %>% pull(gene_name)
})))

study2_gene_pval_summary <- data.frame(
  gene_name = study2_significant_genes,
  highest_pval = sapply(study2_significant_genes, function(gene) {
    pvals <- unlist(lapply(study2_comparisons, function(df) {  # Extract all p-values for the gene from all comparisons
      if (gene %in% df$gene_name) {
        df %>% filter(gene_name == gene) %>% pull(P.Value)
      } else {
        NA  # If the gene is not found, return NA
      }
    }))
    min(pvals, na.rm = TRUE)  # Take the min, ignoring NA values
  })
)

# Sort the genes by their min p-value in ascending order
study2_sorted_genes <- study2_gene_pval_summary %>% arrange(highest_pval)

# Extract logFC values for significant genes across all comparisons
study2_logfc_matrix <- sapply(study2_comparisons, function(df) {
  # Match significant genes and extract logFC
  match_order <- match(study2_sorted_genes$gene_name, df$gene_name)
  df$logFC[match_order]
})

# Assign row names as significant genes
rownames(study2_logfc_matrix) <- study2_sorted_genes$gene_name

study2_logfc_matrix_100 <- head(study2_logfc_matrix, 100)

# Compute the distance matrix
study2_dist_genes <- dist(study2_logfc_matrix, method = "euclidean")
study2_dist_genes_100 <- dist(study2_logfc_matrix_100, method = "euclidean")

# Perform hierarchical clustering using Ward.D2 -clustering is based on min the variation within the clusters
study2_hclust_genes <- hclust(study2_dist_genes, method = "ward.D2")
study2_hclust_genes_100 <- hclust(study2_dist_genes_100, method = "ward.D2")

#jpeg("HierarchialClustering_study2genes_Pval0.001.jpeg", width = 2500, height = 1800, res = 300)
plot(study2_hclust_genes, main = "Hierarchical Clustering of study2 genes at Pval 0.01", xlab = "", sub = "")
#dev.off()

#jpeg("HierarchialClustering_top100genes_Pval0.01.jpeg", width = 4500, height = 1800, res = 300)
plot(study2_hclust_genes_100, main = "Hierarchial Clustering of top 100 genes at Pval 0.01", xlab = "", sub = "")
#dev.off()

S2_breaks <- c(seq(-2, -0.01, length.out = 25), 0, seq(0.01, 2, length.out = 25))

#jpeg("HC_study2genes_Pval0.01.jpeg", width = 3800, height = 3000, res = 300)
pheatmap(
  study2_logfc_matrix,                  # Input matrix
  cluster_rows = study2_hclust_genes,          # Use the dendrogram for rows
  cluster_cols = T,                
  color = colorRampPalette(c("blue", "white", "red"))(50),  # Color gradient
  border_color = NA,
  breaks = S2_breaks,
  main = "Heatmap of genes at Pval 0.01 across study2 comparisons", 
  fontsize_col = 10,                   
  fontsize_row = 6,                    
  labels_col = colnames(study2_logfc_matrix),  
  labels_row = rownames(study2_logfc_matrix), 
  legend = TRUE,                
  angle_col = 45,
  treeheight_row = 200
)
#dev.off()

#jpeg("HC_study2_top100genes_pval0.01.jpeg", width = 3800, height = 2500, res = 300)
pheatmap(
  study2_logfc_matrix_100,                  # Input matrix
  cluster_rows = study2_hclust_genes_100,          # Use the dendrogram for rows
  cluster_cols = FALSE,                
  color = colorRampPalette(c("blue", "white", "red"))(50),  # Color gradient
  border_color = NA,
  main = "Heatmap of top 100 genes at Pval 0.01 across study2 comparisons", 
  fontsize_col = 10,                   
  fontsize_row = 6,                    
  labels_col = colnames(study2_logfc_matrix_100),  
  labels_row = rownames(study2_logfc_matrix_100), 
  legend = TRUE,                
  angle_col = 45,
  treeheight_row = 200
)
#dev.off()

subset <- study2_logfc_matrix[rownames(study2_logfc_matrix) %in% rownames(logfc_matrix), ]

pheatmap(
  subset,                  # Input matrix
  cluster_rows = T,         
  cluster_cols = F,                
  color = colorRampPalette(c("blue", "white", "red"))(50),  # Color gradient
  border_color = NA,
  main = "Heatmap of filtered study2 genes", 
  fontsize_col = 10,                   
  fontsize_row = 6,                    
  labels_col = colnames(subset),  
  labels_row = rownames(subset), 
  legend = TRUE,                
  angle_col = 45
)

##-----------------------------------------------Self-organizing maps - study2 ------------------------------------
library(som)
study2_som = som(study2_logfc_matrix, xdim=4, ydim=3, topol="rect", neigh="bubble", init = "linear")

study2_som_cluster <- study2_som$visual
rownames(study2_som_cluster) <- rownames(study2_logfc_matrix)
#write.xlsx(study2_som_cluster, "SOM_3.xlsx")

#generate plot
#jpeg("SOM_12clusters.jpeg", width = 4500, height = 3500, res = 300)
plot(study2_som, sdbar=1, ylim=c(-3,3), color=T)
#dev.off()

cluster_12 <- study2_som_cluster %>% filter(x == 3 & y == 0)
cluster_8 <- study2_som_cluster %>% filter(x == 3  & y == 1)
cluster_4 <- study2_som_cluster %>% filter(x == 3 & y == 2)
cluster_11 <- study2_som_cluster %>% filter(x == 2 & y == 0)
cluster_7 <- study2_som_cluster %>% filter(x == 2 & y == 1)

combined_cluster_updown <- rbind(cluster_4, cluster_7, cluster_8, cluster_11, cluster_12)


# filter the logfc matrix based on the genes found in the updown pattern cluster
filtered_matrix_updown <- study2_logfc_matrix[rownames(study2_logfc_matrix) %in% rownames(combined_cluster_updown), , drop = FALSE]

#jpeg("Clusters_updown_pattern.jpeg", width = 3500, height = 2500, res = 300)
pheatmap(
  filtered_matrix_updown,                
  cluster_rows = T,          
  cluster_cols = F,                
  color = colorRampPalette(c("blue", "white", "red"))(50),  
  border_color = NA,
  main = "Heatmap of genes at Pval 0.01 with similar expression patterns", 
  fontsize_col = 10,                   
  fontsize_row = 6,                    
  labels_col = colnames(filtered_matrix_updown),  
  labels_row = rownames(filtered_matrix_updown), 
  legend = TRUE,                
  angle_col = 45
)
#dev.off()

cluster_1 <- study2_som_cluster %>% filter(x == 0 & y == 2)
cluster_5 <- study2_som_cluster %>% filter(x == 0 & y == 1)
cluster_2 <- study2_som_cluster %>% filter(x == 1 & y == 2)

combined_cluster_downup <- rbind(cluster_1, cluster_2, cluster_5)

# filter the logfc matrix based on the genes found in the downup pattern cluster
filtered_matrix_downup <- study2_logfc_matrix[rownames(study2_logfc_matrix) %in% rownames(combined_cluster_downup), , drop = F]

#jpeg("Clusters_downup_pattern.jpeg", width = 3500, height = 2500, res = 300)
pheatmap(
  filtered_matrix_downup,                
  cluster_rows = T,          
  cluster_cols = F,                
  color = colorRampPalette(c("blue", "white", "red"))(50),  
  border_color = NA,
  breaks = c(seq(-1, -0.01, length.out = 25), 0, seq(0.01, 1, length.out = 25)),
  main = "Heatmap of genes at Pval 0.01 with similar expression patterns", 
  fontsize_col = 10,                   
  fontsize_row = 6,                    
  labels_col = colnames(filtered_matrix_downup),  
  labels_row = rownames(filtered_matrix_downup), 
  legend = TRUE,                
  angle_col = 45
)
#dev.off()

cluster_3 <- study2_som_cluster %>% filter(x == 2 & y == 2)
cluster_6 <- study2_som_cluster %>% filter(x == 1 & y == 1)
cluster_9 <- study2_som_cluster %>% filter(x == 0 & y == 0)
cluster_10 <- study2_som_cluster %>% filter(x == 1 & y == 0)

cluster_1 <- cluster_1 %>% mutate(Cluster = 1)
cluster_2 <- cluster_2 %>% mutate(Cluster = 2)
cluster_3 <- cluster_3 %>% mutate(Cluster = 3)
cluster_4 <- cluster_4 %>% mutate(Cluster = 4)
cluster_5 <- cluster_5 %>% mutate(Cluster = 5)
cluster_6 <- cluster_6 %>% mutate(Cluster = 6)
cluster_7 <- cluster_7 %>% mutate(Cluster = 7)
cluster_8 <- cluster_8 %>% mutate(Cluster = 8)
cluster_9 <- cluster_9 %>% mutate(Cluster = 9)
cluster_10 <- cluster_10 %>% mutate(Cluster = 10)
cluster_11 <- cluster_11 %>% mutate(Cluster = 11)
cluster_12 <- cluster_12 %>% mutate(Cluster = 12)


combined_df <- bind_rows(cluster_1, cluster_2, cluster_3, cluster_4, cluster_5, cluster_6, cluster_7, cluster_8, cluster_9, cluster_10, cluster_11, cluster_12)
#write.xlsx(combined_df, "SOM_study2.xlsx", rowNames = T)

#Hallmark_Pathway_ud <- read.delim("MSigDB_Hallmark_2020_table (1).txt", sep = '\t')
#Wiki_Pathway_du <- read.delim("WikiPathways_2024_Human_table (1).txt", sep = '\t')
#kegg_Pathway_du <- read.delim("KEGG_2021_Human_table (1).txt", sep = '\t')
