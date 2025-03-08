# New correlation
library(Hmisc)
library(tidyverse)
library(truncnorm)
library(car)
library(ggplot2)
library(ggh4x)
library(grid)

#-------- Generate random data ---------
set.seed(123)  # For reproducibility

n_participants <- 400
n_genes <- 10

sample_id <- sample(paste0("Sample_", sprintf("%03d", 1:n_participants)), n_participants, replace = FALSE)
sex <- sample(c(0, 1), n_participants, replace = TRUE) # 0=male, 1=female
Age <- sample(c(40:60), n_participants, replace = TRUE)
Weight <- sample(c(50:90), n_participants, replace = TRUE)
Ethnicity <- sample(c(1, 2, 3, 4, 5), n_participants, replace = TRUE) #1=A, 2=B, 3=C, 4=D, 5=E
blood_sugar <- sample(c(5:12), n_participants, replace = TRUE)


# Initialize matrix for gene expression
gene_data <- matrix(nrow = n_participants, ncol = n_genes)

# Generate gene expression values based on Sex
for (i in 1:n_participants) {
  if (sex[i] == 0) {
    gene_data[i, ] <- rtruncnorm(n_genes, a = 0, b = 200, mean = 60, sd = 30)
  } else {
    gene_data[i, ] <- rtruncnorm(n_genes, a = 0, b = 200, mean = 40, sd = 20)
  }
}

# Convert to dataframe and add Sex column
colnames(gene_data) <- paste0("Gene_", 1:n_genes)  # Name genes
df <- data.frame(Sampe_ID = sample_id, Sex = sex, Age = Age, 
                 Weight = Weight, Ethnicity = Ethnicity, 
                 Blood_sugar = blood_sugar, gene_data)  # Make the datframe

genes_log10 <- log10(df[,7:ncol(df)])
df_log10 <- cbind(df[1:6], genes_log10)
rownames(df_log10) <- df_log10$Sampe_ID
df_log10 <- df_log10[,-1]

# Correlation analysis for all columns
result_All <- rcorr(as.matrix(df_log10),type = "pearson")
r_values_All <- as.data.frame(as.table(result_All$r))
p_values_All <- as.data.frame(as.table(result_All$P))
n_values_All <- as.data.frame(as.table(result_All$n))

r_p_values_All <- cbind(r_values_All, p_values_All$Freq, n_values_All$Freq)

#Rename the columns of correlation table
colnames(r_p_values_All) <- c("Gene", "Clinical_info", "r_value", "p_value", "n_value")

#Remove the unwanted rows from the correlation table
gene_names <- colnames(df_log10[6:ncol(df_log10)]) #Gene names list


r_p_values_1_All <- subset(r_p_values_All, r_p_values_All$Gene %in% gene_names) 
r_p_values_2_All <- subset(r_p_values_1_All, !(r_p_values_1_All$Clinical_info %in% gene_names)) 
r_p_values_3_All <- r_p_values_2_All %>% mutate("-log10 Pvalue" = -log10(p_value),
                                                ID = "_All")

#For getting samples in each clinical outcome
second_axis_info_All <- r_p_values_3_All[,c("Clinical_info","n_value")] %>% unique() 
colnames(second_axis_info_All)[colnames(second_axis_info_All) == "n_value"] <- "n_value_All"


#Sort the names
r_p_values_3_All$Gene <- factor(r_p_values_3_All$Gene, levels = rev(gene_names))



##Only male samples
df_log10_male <- subset(df_log10, df_log10$Sex == 0)


# Correlation analysis for all columns
result_male <- rcorr(as.matrix(df_log10_male),type = "pearson")
r_values_male <- as.data.frame(as.table(result_male$r))
p_values_male <- as.data.frame(as.table(result_male$P))
n_values_male <- as.data.frame(as.table(result_male$n))

r_p_values_male <- cbind(r_values_male, p_values_male$Freq, n_values_male$Freq)

colnames(r_p_values_male) <- c("Gene", "Clinical_info", "r_value", "p_value", "n_value")


#Remove the unwanted rows from the correlation table
r_p_values_1_male <- subset(r_p_values_male, r_p_values_male$Gene %in% gene_names) 
r_p_values_2_male <- subset(r_p_values_1_male, !(r_p_values_1_male$Clinical_info %in% gene_names)) 
r_p_values_3_male <- r_p_values_2_male %>% mutate("-log10 Pvalue" = -log10(p_value),
                                                ID = "_male")

#For getting samples in each clinical outcome
second_axis_info_male <- r_p_values_3_male[,c("Clinical_info","n_value")] %>% unique() 
colnames(second_axis_info_male)[colnames(second_axis_info_male) == "n_value"] <- "n_value_male"


#Sort the names
r_p_values_3_male$Gene <- factor(r_p_values_3_male$Gene, levels = rev(gene_names))





##Only female samples
df_log10_female <- subset(df_log10, df_log10$Sex == 1)


# Correlation analysis for all columns
result_female <- rcorr(as.matrix(df_log10_female),type = "pearson")
r_values_female <- as.data.frame(as.table(result_female$r))
p_values_female <- as.data.frame(as.table(result_female$P))
n_values_female <- as.data.frame(as.table(result_female$n))

r_p_values_female <- cbind(r_values_female, p_values_female$Freq, n_values_female$Freq)

colnames(r_p_values_female) <- c("Gene", "Clinical_info", "r_value", "p_value", "n_value")


#Remove the unwanted rows from the correlation table
r_p_values_1_female <- subset(r_p_values_female, r_p_values_female$Gene %in% gene_names) 
r_p_values_2_female <- subset(r_p_values_1_female, !(r_p_values_1_female$Clinical_info %in% gene_names)) 
r_p_values_3_female <- r_p_values_2_female %>% mutate("-log10 Pvalue" = -log10(p_value),
                                                  ID = "_female")

#For getting samples in each clinical outcome
second_axis_info_female <- r_p_values_3_female[,c("Clinical_info","n_value")] %>% unique() 
colnames(second_axis_info_female)[colnames(second_axis_info_female) == "n_value"] <- "n_value_female"


#Sort the names
r_p_values_3_female$Gene <- factor(r_p_values_3_female$Gene, levels = rev(gene_names))


# Combine the correlation tables and make the plot
r_p_values_combined <- rbind(r_p_values_3_All, r_p_values_3_male, r_p_values_3_female)

r_p_values_combined <- r_p_values_combined%>% mutate(Clinical_info_2 = paste0(Clinical_info, ID),
                                                     Gene_2 = paste0(Gene, ID))


r_p_values_combined <- r_p_values_combined %>% arrange(desc(Gene), Clinical_info, desc(Gene_2))


List_clinical <- c("Sex", "Age", "Ethnicity", "Weight", "Blood_sugar")


r_p_values_combined$Clinical_info <- factor(r_p_values_combined$Clinical_info, 
                                            levels = rev(List_clinical))
r_p_values_combined$Gene_2 <- factor(r_p_values_combined$Gene_2, levels = unique(r_p_values_combined$Gene_2))
r_p_values_combined$Gene <- factor(r_p_values_combined$Gene, levels = unique(r_p_values_combined$Gene))



second_axis_info_list <- list(second_axis_info_All, second_axis_info_male, second_axis_info_female)
second_axis_info_data <- second_axis_info_list %>% purrr::reduce(merge, by = "Clinical_info")
second_y_axis_info_data <- second_axis_info_data %>% mutate(all_n = paste0("(",n_value_All,") ",
                                                                           "(",n_value_male,") ",
                                                                           "(",n_value_female,")"))






r_p_values_combined_filtered <- r_p_values_combined %>% mutate(r_value = ifelse((r_value <= 0.1 & 
                                                                                   r_value >= -0.1), 0, r_value),
                                                               p_value = ifelse(r_value >= 0.1, 1, p_value))




p <- ggplot(data = r_p_values_combined_filtered, 
            aes(x = Gene_2, fill = r_value, 
                y = Clinical_info, size = `-log10 Pvalue`)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  scale_fill_gradientn(colours = c("blue","white", "red"), 
                       limit = c(-0.45, 0.2), 
                       breaks = c(-0.45, -0.3, 0, 0.1, 0.2), 
                       labels=c("-0.5", "-0.3", "0", "0.1", "0.2")) + #Set the color change according to your correlation
  scale_size_continuous(breaks = c(0, 0.5, 1.5, 2, 5), 
                        labels = c("<1","1","2.5","5","10"), 
                        range = c(0,5)) + #Set the size, change according to your log10 pvalues
  facet_wrap( ~ Gene, scales = "free_x", 
              nrow = 1, strip.position = "bottom") + #Group the correlations for each gene
  theme(panel.background = element_rect(fill = NA, color = "black", linewidth = 0.2), 
        plot.background = element_rect(fill = NA, color = NA),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.2),
        panel.grid.major.x = element_line(linewidth = 0.1, colour = "grey80"),
        panel.grid.major.y = element_line(linewidth = 0.1, colour = "grey80"),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_text(angle = 0, size = 8, vjust = 0.4,hjust = 1, color = "black"),
        axis.text.y.right = element_text(angle = 0, size = 7, vjust = 0.4,hjust = 1, color = "black"),
        axis.title.y.right = element_text(angle = -90, size = 7, 
                                          vjust = 0.5,hjust = 0.5, color = "black", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 10)), 
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank(),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.2),
        axis.ticks.x = element_blank(),
        legend.title = element_text(angle = 0, vjust = 0.5, hjust = 0.5,size = 8, colour = "black"),
        legend.key = element_rect(colour = "white"),
        strip.text = element_text(size = 12, color = "black", angle = 90, hjust = 1, 
                                  margin = margin(t = 14, r = 0, b = 10, l = 0)), 
        strip.background = element_rect(fill = NA),
        panel.spacing = unit(0.5, "mm")
  ) +
  guides(
         y.sec = guide_axis_manual(breaks = second_y_axis_info_data$Clinical_info, 
                                   title = "# of Samples in the correlation \n (Male) (Female) (All)",
                                   labels = second_y_axis_info_data$all_n, label_size = 7),
         fill = guide_colorbar(title = "Correlation \n(pearson)", title.position = "top", 
                               order = 1, barwidth = 0.6, barheight = 4),
         size = guide_legend(title = "-log10 \np-value", 
                             title.position = "top",  order = 2, )) + #lables
  annotation_custom(grob = rectGrob(gp = gpar(fill = "skyblue3", col = "black")), 
                    xmin = 1 - 0.5, xmax = 1 + 0.5, ymin = 0.15, ymax = 0.3) +
  annotation_custom(grob = rectGrob(gp = gpar(fill = "orange1", col = "black")), 
                    xmin = 2 - 0.5, xmax = 2 + 0.5, ymin = 0.15, ymax = 0.3) +
  annotation_custom(grob = rectGrob(gp = gpar(fill = "darkseagreen", col = "black")), 
                    xmin = 3 - 0.5, xmax = 3 + 0.5, ymin = 0.15, ymax = 0.3) + # Color annotation at the bottom
  coord_cartesian(clip = "off") #To remove background colors


ggsave("Correlation_balloon_plot.pdf", plot = p, device = "pdf", width = 18, height = 9, units = "cm", ,bg = "transparent")

