##### Libraries #####

if(!require(data.table)) install.packages("data.table", repos = "http://cran.us.r-project.org")
if(!require(cowplot)) install.packages("cowplot", repos = "http://cran.us.r-project.org")
if(!require(vegan)) install.packages("vegan", repos = "http://cran.us.r-project.org")
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(kableExtra)) install.packages("kableExtra ", repos = "http://cran.us.r-project.org")
if(!require(reshape2)) install.packages("reshape2 ", repos = "http://cran.us.r-project.org")

library(data.table)
library(cowplot)
library(vegan)
library(tidyverse)
library(caret)
library(kableExtra)
library(reshape2)


##### Read file #####

# the file were downloaded from here:
# https://www.kaggle.com/datasets/antaresnyc/human-gut-microbiome-with-asd?select=GSE113690_Autism_16S_rRNA_OTU_assignment_and_abundance.csv

Raw <- fread("GSE113690_Autism_16S_rRNA_OTU_assignment_and_abundance.csv", 
             stringsAsFactors = T)


##### Data organization ######

# Create Taxonomy data.frame, with OTU and Taxonomy columns
Taxonomy <- Raw %>% 
  select(OTU, taxonomy)

# Remove Taxonomy and OTU columns, 
# transpose the results and 
# rename the columns with the OTU column to cross-reference the taxonomy table.
# The rownames are then transformed into a column and the tags are 
# recoded as ASD or TD (According to the supplementary material of the paper).
# Any bacteria without reads are then removed from the object.
Data <- Raw %>% 
  select(-taxonomy, -OTU) %>% 
  t() %>% 
  data.frame() %>% 
  `colnames<-`(Raw$OTU) %>% 
  rownames_to_column("target") %>% 
  mutate(target = str_remove_all(target, "\\d"),
         target = factor(target),
         target = fct_recode(target, ASD = "B", TD = "A"),
  ) %>%
  select(target, where(~ is.numeric(.x) && sum(.x) != 0))

# A new data.frame is created with the target column and
# this column is removed from the Data object.
Target <- select(Data, target)
Data <- select(Data, -target)

# The bacteria present in the taxonomy object is filtered based on the remaining 
# bacteria in the dataset.
# The taxonomy column is separated in 8 columns, with each column representing 
# a taxonomy rank.
Taxonomy <- Taxonomy %>% 
  filter(OTU %in% colnames(Data)) %>% 
  separate(taxonomy, sep = ";_", 
           into = c("domain",
                    "kingdom",
                    "phylum", 
                    "class", 
                    "order",
                    "family", 
                    "genus", 
                    "species"))




######## Figure 1 - Number of samples in each group ########

Target %>% 
  count(target) %>% 
  ggplot(aes(x = target, y = n, fill = target, label = n)) +
  geom_col() + 
  geom_text(vjust = 1) +
  theme_bw() +
  xlab("Condition") +
  ylab("Number of samples") +
  theme(legend.position = "none")





######## Figure 2 - Stacked bar-plot ########

# All three objects are united into one.
# first the data is converted to relative abundance,
# then it is melted to long format, suitable to ploting,
# then its joined to the taxonomy table
# and finally with the target table.
# We can summarise the resulting object to create compositional bar plots of different taxa ranks.

FullData <- Data %>% 
  apply(MARGIN = 1, function(x) { (x / sum(x))*100 }) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column("Sample") %>% 
  melt(id.vars = "Sample", variable.name = "OTU", value.name = "Abundance") %>% 
  left_join(Taxonomy, by = "OTU") %>% 
  left_join(rownames_to_column(Target, "Sample"), by = "Sample") 

FullData %>% 
  group_by(Sample, target, phylum) %>% 
  summarise(Abundance = sum(Abundance), .groups = "drop_last") %>% 
  ggplot(aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_col() +
  facet_grid(~ target, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(2, 'mm'),
        axis.text.x = element_blank()) +
  guides(fill=guide_legend(ncol=3)) +
  ylab("Relative Abundance (%)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0)))




######## Figure 3 - PCA ########

Data.pca <- prcomp(Data)

PCA_1_2 <- Data.pca$x %>% 
  data.frame() %>% 
  mutate(target = Target$target) %>% 
  ggplot(aes(x = PC1, y = PC2, color = target)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab(paste0("PC 1 (", round(100 * Data.pca$sdev[1]^2/sum(Data.pca$sdev^2), 1), "%)")) +
  ylab(paste0("PC 2 (", round(100 * Data.pca$sdev[2]^2/sum(Data.pca$sdev^2), 1), "%)"))

PCA_3_4 <- Data.pca$x %>% 
  data.frame() %>% 
  mutate(target = Target$target) %>% 
  ggplot(aes(x = PC3, y = PC4, color = target)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab(paste0("PC 3 (", round(100 * Data.pca$sdev[3]^2/sum(Data.pca$sdev^2), 1), "%)")) +
  ylab(paste0("PC 4 (", round(100 * Data.pca$sdev[4]^2/sum(Data.pca$sdev^2), 1), "%)"))

# Plot two graphs in a grid, without legend
# then put these two in a new grid, over the legend
plot_grid(PCA_1_2 + theme(legend.position = "none"), 
          PCA_3_4 + theme(legend.position = "none"), 
          labels = c("A", "B"),
          ncol = 2) %>% 
  plot_grid(get_legend(PCA_1_2), 
            ncol = 1,
            rel_heights = c(0.9, 0.1))


######## Figure 4 - PCoA ########


# function to perform PCoA, from the package vegan
# to calculate de jaccard dissimilarity, it is important to set the binary option to TRUE
# the option metaMDSdist automatic transforms the data if necessary and decreases the bias
# caused by the sparse data.
Data.pco.jaccard <- capscale(Data ~ 1, distance = "jaccard", binary = T, metaMDSdist = TRUE, add = TRUE)

# function scores extracts the coordinates from the ordination result
# Plot axis 1 and 2
PCO.jaccard_1_2 <- scores(Data.pco.jaccard, display = "sites") %>% 
  as.data.frame() %>% 
  mutate(target = Target$target) %>% 
  ggplot(aes(x = MDS1, y = MDS2, color = target)) +
  geom_point() +
  theme_bw() +
  stat_ellipse() +
  xlab(paste0("PCoA 1 (", round(100 * Data.pco.jaccard$CA$eig[1]/sum(Data.pco.jaccard$CA$eig), 1), "%)")) +
  ylab(paste0("PCoA 2 (", round(100 * Data.pco.jaccard$CA$eig[2]/sum(Data.pco.jaccard$CA$eig), 1), "%)")) + 
  theme(legend.position = "bottom")

# Plot axis 3 and 4
PCO.jaccard_3_4 <- scores(Data.pco.jaccard, display = "sites", choices = c(3, 4)) %>% 
  as.data.frame() %>% 
  mutate(target = Target$target) %>% 
  ggplot(aes(x = MDS3, y = MDS4, color = target)) +
  geom_point() +
  theme_bw() +
  stat_ellipse() +
  xlab(paste0("PCoA 3 (", round(100 * Data.pco.jaccard$CA$eig[3]/sum(Data.pco.jaccard$CA$eig), 1), "%)")) +
  ylab(paste0("PCoA 4 (", round(100 * Data.pco.jaccard$CA$eig[4]/sum(Data.pco.jaccard$CA$eig), 1), "%)")) 


# PCoA based on Bray-Curtis dissimilarity
Data.pco.bray <- capscale(Data ~ 1, distance = "bray", binary = F, metaMDSdist = TRUE, add = TRUE)

# Plot axis 1 and 2
PCO.bray_1_2 <- scores(Data.pco.bray, display = "sites") %>% 
  as.data.frame() %>% 
  mutate(target = Target$target) %>% 
  ggplot(aes(x = MDS1, y = MDS2, color = target)) +
  geom_point() +
  theme_bw() +
  stat_ellipse() +
  xlab(paste0("PCoA 1 (", round(100 * Data.pco.bray$CA$eig[1]/sum(Data.pco.bray$CA$eig), 1), "%)")) +
  ylab(paste0("PCoA 2 (", round(100 * Data.pco.bray$CA$eig[2]/sum(Data.pco.bray$CA$eig), 1), "%)"))

# Plot axis 3 and 4
PCO.bray_3_4 <- scores(Data.pco.bray, display = "sites", choices = c(3, 4)) %>% 
  as.data.frame() %>% 
  mutate(target = Target$target) %>% 
  ggplot(aes(x = MDS3, y = MDS4, color = target)) +
  geom_point() +
  theme_bw() +
  stat_ellipse() +
  xlab(paste0("PCoA 3 (", round(100 * Data.pco.bray$CA$eig[3]/sum(Data.pco.bray$CA$eig), 1), "%)")) +
  ylab(paste0("PCoA 4 (", round(100 * Data.pco.bray$CA$eig[4]/sum(Data.pco.bray$CA$eig), 1), "%)")) 


# Full plot with all 4 graphs:
plot_grid(PCO.jaccard_1_2 + theme(legend.position = "none"), 
          PCO.jaccard_3_4 + theme(legend.position = "none"), 
          PCO.bray_1_2 + theme(legend.position = "none"), 
          PCO.bray_3_4 + theme(legend.position = "none"), 
          ncol = 2, 
          labels = c("A", "B", "C", "D")) %>% 
  plot_grid(get_legend(PCO.jaccard_1_2), 
            ncol = 1,
            rel_heights = c(0.95, 0.05))






######## Modeling approach ########

set.seed(2, sample.kind = "Rounding") # if using R 3.5 or earlier, remove the sample.kind argument

# Create index for each dataset
test_index <- createDataPartition(Target$target, times = 1, p = 0.2, list = FALSE)

# Create train set
train_set <- Data[-test_index, ]
train_set_y <- Target$target[-test_index] 

# Create test set
test_set <- Data[test_index, ]
test_set_y <- Target$target[test_index] 


# Remove bacteria not present in the train set
train_set <- train_set[,colSums(train_set) != 0]
test_set <- test_set[, colnames(test_set) %in% colnames(train_set)]

# Transforme to presence/ausence
train_set_PA <- decostand(train_set, "pa")
test_set_PA <- decostand(test_set, "pa")

# Train method for the models: 5-fold cross validation
control <- trainControl(method = "cv", number = 5, p = 0.8)



######## Model building ########

Models <- c("Rborist", "xgbTree", "pcaNNet", "glm", "naive_bayes", "knn")

# These objects will save the results in the loop
fits <- NULL
fits_PA <- NULL
Y_hats <- NULL
Y_hats_PA <- NULL
Performance <- NULL
for(i in Models) {
  
  #### Train models in the raw data
  
  # Set seed
  # set.seed(2) #if you are using R 3.5 or earlier
  set.seed(2, sample.kind = "Rounding") 
  
  # keep track of starts and end times
  start_time <- Sys.time()
  fits[[i]] <- train(x = train_set,
                     y = train_set_y,
                     method = i, 
                     trControl = control)
  end_time <- Sys.time()
  
  # make predictions in the test set
  Y_hats[[i]] <- predict(fits[[i]], test_set)
  
  # Calculate and save performancy metrics
  CM <- confusionMatrix(Y_hats[[i]],
                        test_set_y, 
                        positive = "ASD")
  
  
  
  #### Train models in the PA data
  
  # Set seed
  # set.seed(2) #if you are using R 3.5 or earlier
  set.seed(2, sample.kind = "Rounding") 
  
  # keep track of starts and end times
  start_time_PA <- Sys.time()
  fits_PA[[i]] <- train(x = train_set_PA,
                        y = train_set_y,
                        method = i, 
                        trControl = control)
  end_time_PA <- Sys.time()
  
  # make predictions in the test set
  Y_hats_PA[[i]] <- predict(fits_PA[[i]], test_set_PA)
  
  # Calculate and save performancy metrics
  CM_PA <- confusionMatrix(Y_hats_PA[[i]],
                           test_set_y, 
                           positive = "ASD")
  
  
  
  # final performance data:
  Performance[[i]] <- list(Accuracy = CM$overall[["Accuracy"]],
                           Accuracy_PA = CM_PA$overall[["Accuracy"]],
                           Sensitivity = CM$byClass[["Sensitivity"]],
                           Sensitivity_PA = CM_PA$byClass[["Sensitivity"]],
                           Specificity = CM$byClass[["Specificity"]],
                           Specificity_PA = CM_PA$byClass[["Specificity"]],
                           Time = as.numeric(difftime(end_time, start_time, units = "secs")),
                           Time_PA = as.numeric(difftime(end_time_PA, start_time_PA, units = "secs")))
}



######## Results ########

# Performance table:

# transform the Performance list into a data frame,
# round all numbers to 3 decimal places
# use kable to create a table
rbindlist(Performance, idcol = "model")%>% 
  mutate(across(where(is.numeric), ~ round(.x, 3) )) %>%  
  kbl(booktabs = T, 
      linesep = "", 
      align = "c",
      col.names = c("", "Raw", "P/A", "Raw", "P/A", "Raw", "P/A", "Raw", "P/A"),
      caption = "Performance metrics of the six different models trained in the microbiome dataset.") %>% 
  add_header_above(c("", "Accuracy" = 2, "Sensitivity" = 2, "Specifcity" = 2, "Time (seconds)" = 2)) %>% 
  kable_styling(latex_options = c("striped", "hold_position"))


#Figure 5 - Performance:

# transform the Performance list into a data frame,
# convert to long format
# then organize the columns to plot
rbindlist(Performance, idcol = "model") %>% 
  melt(id.vars = "model") %>% 
  mutate(Data = case_when(str_detect(variable, "_PA") ~ "PA",
                          TRUE ~ "Raw"),
         variable = str_remove_all(variable, "_PA"),
         variable = case_when(str_detect(variable, "Time") ~ paste(variable, "(seconds)"),
                              TRUE ~ paste(variable, "(%)")),
         model = factor(model, levels = unique(model))) %>% 
  ggplot(aes(x = model, y = value, color = Data, group = Data)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ variable, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom")



######## Conclusion ########

# confusion matrix for knn in the PA dataset:
confusionMatrix(Y_hats_PA$knn,
                test_set_y, 
                positive = "ASD")$table