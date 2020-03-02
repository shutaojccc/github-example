#' ---
#' title: "Boutros Lab Onboarding R Training"
#' author: "Shu Tao"
#' date: "2/22/2020"
#' output: word_document
#' ---
#' 
### RTraining.R ####################################################################################
# This R script was used to complete R training in the '1. Two Sample Tests' and '2. Plotting' folder.

### PREAMBLE #######################################################################################
# Load Boutros Lab R packages.
library(BoutrosLab.statistics.general);
library(BoutrosLab.plotting.general);

# General parameters.
date <- Sys.Date();

# Q1 Answers:
# Set working directory.
setwd('/Users/sthoya/Downloads/RTraining/1. Two Sample Tests/');

# Read the input file: AHR-test-file.txt.
AHR.test.input <- read.table(
  file = 'AHR-test-file.txt',
  header = TRUE,
  sep = '\t',
  row.names = 1 # Set the first column of the input file as row names.
);

# Perform a t-test between control and treated.
control.data <- as.numeric(AHR.test.input$Control);
treated.data <- as.numeric(AHR.test.input$Treated);
t.test(control.data,
       treated.data,
       alternative = 'two.sided',
       var.equal = FALSE,
       paired = TRUE #Paired t-test because each sample was measured in both control and treatment conditions.
       );

# Perform a wilcoxon test between control and treated.
wilcox.test(control.data,
            treated.data,
            alternative = 'two.sided',
            paired = TRUE
            );

# Calculate a fold-change between control and treated.
AHR.test.input$foldchange <- treated.data/control.data;
AHR.test.input; # Display the fold change result.


# Q2 Answers:
# Read the input data for each tumor subtype.
# Tumor subtype A.
tumor.a.expression <- read.table(
  file = 'input1.txt',
  header = TRUE,
  sep = '\t',
  row.names = 1
  );

# Tumor subtype B.
tumor.b.expression <- read.table(
  file = 'input2.txt',
  header = TRUE,
  sep = '\t',
  row.names = 1
  );

# Combine data from two tumor subtypes.
# Method 1: sort each tumor data individually and use cbind function to combine them.
tumor.a.expression.sort <- tumor.a.expression[order(row.names(tumor.a.expression),
                                                    decreasing = FALSE
                                                    ),
                                              ];

tumor.b.expression.sort <- tumor.b.expression[order(row.names(tumor.b.expression),
                                                    decreasing = FALSE
                                                    ),
                                              ];

tumor.ab.expression.sort.cbind <- cbind(tumor.a.expression.sort,
                                        tumor.b.expression.sort
                                        );

# Method 2: use only the merge function.
tumor.ab.expression.merge <- merge(tumor.a.expression,
                                   tumor.b.expression,
                                   by = 'row.names' # Merge by row names (i.e., gene IDs).
                                   );

tumor.ab.expression.merge.clean <- data.frame(tumor.ab.expression.merge,
                                              row.names = 1
                                              ); # Use row.names = 1 to assign the 'Row.names' column as row names.

### Function all.equal #############################################################################
# Description: use all.equal() function to verify that both merging methods produce the same result.
# Input variable:
#   tumor.ab.expression.sort.cbind: data frame merged using cbind() function.
#   tumor.ab.expression.merge.clean: data frame merged using merge() function.
# Output variable:
#   compare.merge.result: TRUE indicating same merging results; otherwise difference is reported. 
# After verification, use tumor.ab.expression.sort.cbind in subsequent analysis.
compare.merge.result <- all.equal(tumor.ab.expression.sort.cbind,
                                  tumor.ab.expression.merge.clean
                                  );
compare.merge.result; # Return TRUE.

# Perform a T-Test comparing the first three tumours to the last nine tumours for each gene using a for-loop.
tumor.a.patients <- colnames(tumor.a.expression); # Get patient names (1 through 3) from Tumor subtype A.
tumor.b.patients <- colnames(tumor.b.expression); # Get patient names (4 through 12) from Tumor subtype B.

# Use for-loop to perform T-Test for each gene between two tumor subtypes.
# Write T-Test P-Values to a new column with name 'ttest.p.value'.
for (gene.counter in 1:nrow(tumor.ab.expression.sort.cbind)){
  tumor.a.gene <- as.numeric(tumor.ab.expression.sort.cbind[gene.counter, tumor.a.patients]);
  tumor.b.gene <- as.numeric(tumor.ab.expression.sort.cbind[gene.counter, tumor.b.patients]);
  gene.ttest <- t.test(tumor.a.gene,
                       tumor.b.gene,
                       alternative = 'two.sided',
                       var.equal = FALSE,
                       paired = FALSE
                       );
  tumor.ab.expression.sort.cbind[gene.counter, 'ttest.p.value'] <- gene.ttest$p.value;
  }

# Plot a histogram of the t-test p-values.
hist(tumor.ab.expression.sort.cbind$ttest.p.value,
     main = 'Histogram of T Test P-Values', # Sets histogram main title.
     freq = TRUE, # Defaults to TRUE if and only if breaks of x axis are equidistant.
     xlab = 'P-Value',
     ylab = 'Counts', # Change default ylab name from 'Frequency' to 'Counts' because 'freq = TRUE' displays the counts component of the result.
     ylim = c(0, 80), # Sets Y axis label limits.
     las = 1 # Make y axis labels horizontal.
     );

# Plot a histogram using log10-based p-values.
hist(log10(tumor.ab.expression.sort.cbind$ttest.p.value),
     main = 'Histogram of T Test P-Values',
     freq = TRUE,
     xlab = 'Log10-based P-Value',
     ylab = 'Counts',
     ylim = c(0, 360),
     las = 1
     ); # The P-Values distribution shows that the number of genes with a significant difference (P-Value < 0.05) of average expression between Tumor subtype A and B is very few.


# Q3 Answers:
# Use for-loop to perform Wilcoxon Rank Sum Test for each gene between two tumor subtypes.
# Write Wilcoxon Rank Sum Test P-Values to a new column with name 'wilcoxon.p.value'.
# Write gene expression fold change to a new column with name 'FoldChange'.
for (gene.counter in 1:nrow(tumor.ab.expression.sort.cbind)){
  tumor.a.gene <- as.numeric(tumor.ab.expression.sort.cbind[gene.counter, tumor.a.patients]);
  tumor.b.gene <- as.numeric(tumor.ab.expression.sort.cbind[gene.counter, tumor.b.patients]);
  gene.wilcoxontest <- wilcox.test(tumor.a.gene,
                                   tumor.b.gene,
                                   alternative = 'two.sided',
                                   var.equal = FALSE,
                                   paired = FALSE
                                   );
  tumor.ab.expression.sort.cbind[gene.counter, 'wilcoxon.p.value'] <- gene.wilcoxontest$p.value;
  tumor.ab.expression.sort.cbind[gene.counter, 'FoldChange'] <- median(tumor.b.gene)/median(tumor.a.gene);
  }

# Plot a histogram of the Wilcoxon Rank Sum Test p-values.
hist(tumor.ab.expression.sort.cbind$wilcoxon.p.value,
     main = 'Histogram of Wilcoxon Rank Sum Test P-Values',
     freq = TRUE,
     xlab = 'P-Value',
     ylab = 'Counts',
     ylim = c(0, 80),
     las = 1
     );

#> sum((tumor.ab.expression.sort.cbind$wilcoxon.p.value > 0) & (tumor.ab.expression.sort.cbind$wilcoxon.p.value <= 0.1))
#[1] 62
#> sum((tumor.ab.expression.sort.cbind$wilcoxon.p.value > 0.1) & (tumor.ab.expression.sort.cbind$wilcoxon.p.value <= 0.2))
#[1] 34
#> sum((tumor.ab.expression.sort.cbind$wilcoxon.p.value > 0.2) & (tumor.ab.expression.sort.cbind$wilcoxon.p.value <= 0.3))
#[1] 58
#> sum((tumor.ab.expression.sort.cbind$wilcoxon.p.value > 0.3) & (tumor.ab.expression.sort.cbind$wilcoxon.p.value <= 0.4))
#[1] 46
#> sum((tumor.ab.expression.sort.cbind$wilcoxon.p.value > 0.4) & (tumor.ab.expression.sort.cbind$wilcoxon.p.value <= 0.5))
#[1] 38
#
#> sum((tumor.ab.expression.sort.cbind$wilcoxon.p.value > 0.5) & (tumor.ab.expression.sort.cbind$wilcoxon.p.value <= 0.6))
#[1] 59
#> sum((0.6 == tumor.ab.expression.sort.cbind$wilcoxon.p.value))
#[1] 55
#> sum((tumor.ab.expression.sort.cbind$wilcoxon.p.value > 0.5) & (tumor.ab.expression.sort.cbind$wilcoxon.p.value < 0.6))
#[1] 4
#
#> sum((tumor.ab.expression.sort.cbind$wilcoxon.p.value > 0.6) & (tumor.ab.expression.sort.cbind$wilcoxon.p.value <= 0.7))
#[1] 0
#> sum((tumor.ab.expression.sort.cbind$wilcoxon.p.value > 0.7) & (tumor.ab.expression.sort.cbind$wilcoxon.p.value <= 0.8))
#[1] 66
#> sum((tumor.ab.expression.sort.cbind$wilcoxon.p.value > 0.8) & (tumor.ab.expression.sort.cbind$wilcoxon.p.value <= 0.9))
#[1] 70
#> sum((tumor.ab.expression.sort.cbind$wilcoxon.p.value > 0.9) & (tumor.ab.expression.sort.cbind$wilcoxon.p.value <= 1))
#[1] 67

# Plot a histogram of the Wilcoxon Rank Sum Test using log10-based p-values.
hist(log10(tumor.ab.expression.sort.cbind$wilcoxon.p.value),
     main = 'Histogram of Wilcoxon Rank Sum Test P-Values',
     freq = TRUE,
     xlab = 'Log10-based P-Value',
     ylab = 'Counts',
     ylim = c(0, 360),
     las = 1
     );

# Create a data frame to record P-Values from each test.
t.and.wilcox.p.values <- data.frame(
  P.Value = c(tumor.ab.expression.sort.cbind$ttest.p.value,
              tumor.ab.expression.sort.cbind$wilcoxon.p.value
              ),
  Test = c(rep('T',
               length(tumor.ab.expression.sort.cbind$ttest.p.value)
               ),
           rep('Wilcoxon Rank Sum',
               length(tumor.ab.expression.sort.cbind$wilcoxon.p.value)
               )
           )
  );

# Compare P-Values from T-Test and Wilcoxon Rank Sum Test using boxplot.
plot(P.Value ~ Test,
     data = t.and.wilcox.p.values,
     lwd = 1, # Sets line width.
     main = 'Comparison of P-Values',
     xlab = 'Test',
     ylab = 'P-Value'
     );

### Function stripchart ############################################################################
# Description: use stripchart() function to add scattered points to the boxplot.
# Input variable:
#   t.and.wilcox.p.values: a data frame of P-Values and tests' names.
# Output variable:
#   Added scattered points in the boxplot.
stripchart(P.Value ~ Test,
           data = t.and.wilcox.p.values,
           add = TRUE, # Add points to the current plot.
           method = "jitter", # Jitter points to avoid overplotting.
           vertical = TRUE, # Points are drawn vertically.
           pch = 19,
           cex = 0.2,
           col = 'blue'
           ); # The scattered points added to the box plot show that both tests have similar range (0 to 1) of P-Values. However, T Test produces P-Values that are more evenly scattered, while Wilcoxon Rank Sum Test generates very discrete P-Values.

# Use 'apply' function to generate a vector of T Test P-Values.
### Function calc.ttest.pval #######################################################################
# Description: use calc.ttest.pval() function to calculate T-Test P-Value for each gene.
# Input variable:
#   gene.expr: A numeric vector to represent gene expression in all 12 samples.
# Output variable:
#   gene.ttest.in.function$p.value: T Test P-Value for the gene between two tumor subtypes.
calc.ttest.pval <- function(gene.expr) {
  tumor.a.gene <- as.numeric(gene.expr[tumor.a.patients]);
  tumor.b.gene <- as.numeric(gene.expr[tumor.b.patients]);
  gene.ttest.in.function <- t.test(tumor.a.gene,
                                   tumor.b.gene,
                                   alternative = 'two.sided',
                                   var.equal = FALSE,
                                   paired = FALSE
                                   );
  return(gene.ttest.in.function$p.value);
  }

# Use apply() function to invoke calc.ttest.pval() function to calculate T Test P-Value
# for each row (gene) of the input data frame tumor.ab.expression.sort.cbind.
ttest.p.value.by.apply <- apply(tumor.ab.expression.sort.cbind,
                                1, # Sets 1 to represent row (gene) level.
                                calc.ttest.pval # Invoke calc.ttest.pval() function.
                                );

# Confirm T Test P-Values calclated by calc.ttest.pval() function using 'apply' wrapper
# is the same as those calculated by for-loop earlier.
all.equal(tumor.ab.expression.sort.cbind$ttest.p.value, # T Test P-Values calculated by for-loop earlier.
          unname(ttest.p.value.by.apply) # T Test P-Values calclated by calc.ttest.pval() function here.
          ); # Return TRUE.

# Use 'apply' function to generate a vector of Wilcoxon Rank Sum Test P-Values.
### Function calc.wilcoxontest.pval ################################################################
# Description: use calc.wilcoxontest.pval() function to calculate Wilcoxon Rank Sum Test P-Value for each gene.
# Input variable:
#   gene.expr: A numeric vector to represent gene expression in all 12 samples.
# Output variable:
#   gene.wilcoxontest.in.function$p.value: Wilcoxon Rank Sum Test P-Value for the gene between two tumor subtypes.
calc.wilcoxontest.pval <- function(gene.expr) {
  tumor.a.gene <- as.numeric(gene.expr[tumor.a.patients]);
  tumor.b.gene <- as.numeric(gene.expr[tumor.b.patients]);
  gene.wilcoxontest.in.function <- wilcox.test(tumor.a.gene,
                                               tumor.b.gene,
                                               alternative = 'two.sided',
                                               var.equal = FALSE,
                                               paired = FALSE
                                               );
  return(gene.wilcoxontest.in.function$p.value);
  }

# Use apply() function to invoke calc.wilcoxontest.pval() function to calculate Wilcoxon Rank Sum Test
# P-Value for each row (gene) of the input data frame tumor.ab.expression.sort.cbind.
wilcoxontest.p.value.by.apply <- apply(tumor.ab.expression.sort.cbind,
                                       1, # Sets 1 to represent row (gene) level.
                                       calc.wilcoxontest.pval # Invoke calc.wilcoxontest.pval() function.
                                       );

# Confirm Wilcoxon Rank Sum Test P-Values calclated by calc.wilcoxontest.pval() function
# using 'apply' wrapper is the same as those calculated by for-loop earlier.
all.equal(tumor.ab.expression.sort.cbind$wilcoxon.p.value, # Wilcoxon Rank Sum Test P-Values calculated by for-loop earlier.
          unname(wilcoxontest.p.value.by.apply) # Wilcoxon Rank Sum Test P-Values calclated by calc.wilcoxontest.pval() function here.
          ); # Return TRUE.

# Use 'apply' function to generate a vector of fold change.
### Function calc.foldchange #######################################################################
# Description: use calc.foldchange() function to calculate fold change for each gene.
# Input variable:
#   gene.expr: A numeric vector to represent gene expression in all 12 samples.
# Output variable:
#   foldchange.in.function: fold change for the gene between two tumor subtypes.
calc.foldchange <- function(gene.expr) {
  tumor.a.gene <- as.numeric(gene.expr[tumor.a.patients]);
  tumor.b.gene <- as.numeric(gene.expr[tumor.b.patients]);
  foldchange.in.function <- median(tumor.b.gene)/median(tumor.a.gene);
  return(foldchange.in.function);
}

# Use apply() function to invoke calc.foldchange() function to calculate fold change
# for each row (gene) of the input data frame tumor.ab.expression.sort.cbind.
foldchange.by.apply <- apply(tumor.ab.expression.sort.cbind,
                             1, # Sets 1 to represent row (gene) level.
                             calc.foldchange # Invoke calc.foldchange() function.
                             );

# Confirm fold change calclated by calc.foldchange() function
# using 'apply' wrapper is the same as those calculated by for-loop earlier.
all.equal(tumor.ab.expression.sort.cbind$FoldChange, # Fold change calculated by for-loop earlier.
          unname(foldchange.by.apply) # Fold change calclated by calc.foldchange() function here.
          ); # Return TRUE.


# Q4 Answers:
# Adjust T Test P-Values for multiple-testing using FDR and Bonferonni.
tumor.ab.expression.sort.cbind$ttest.p.value.FDR <- p.adjust(
  tumor.ab.expression.sort.cbind$ttest.p.value,
  method = 'fdr'
  );
tumor.ab.expression.sort.cbind$ttest.p.value.bonferroni <- p.adjust(
  tumor.ab.expression.sort.cbind$ttest.p.value,
  method = 'bonferroni'
  );

par(mfrow=c(1,3)); # Make three P-Values next to one another in a row for comparison.

hist(tumor.ab.expression.sort.cbind$ttest.p.value,
     main = 'Histogram of T Test\nP-Values',
     freq = TRUE,
     xlab = 'P-Value',
     ylab ='Counts',
     ylim = c(0, 500),
     las = 1
     );

hist(tumor.ab.expression.sort.cbind$ttest.p.value.FDR,
     main = 'Histogram of T Test\nP-Values Adjusted by\nFDR',
     freq = TRUE,
     xlab = 'P-Value (FDR)',
     ylab = 'Counts',
     ylim = c(0, 500),
     las = 1
     );

hist(tumor.ab.expression.sort.cbind$ttest.p.value.bonferroni,
     main = 'Histogram of T Test\nP-Values Adjusted by\nBonferroni',
     freq = TRUE,
     xlab = 'P-Value (Bonferroni)',
     ylab = 'Counts',
     ylim = c(0, 500),
     las = 1
     );

par(mfrow=c(1,1)); # Resets the mfrow parameter. 

# Adjust Wilcoxon Rank Sum Test P-Values for multiple-testing using FDR and Bonferonni.
tumor.ab.expression.sort.cbind$wilcoxon.p.value.FDR <- p.adjust(
  tumor.ab.expression.sort.cbind$wilcoxon.p.value,
  method = 'fdr'
  );
tumor.ab.expression.sort.cbind$wilcoxon.p.value.bonferroni <- p.adjust(
  tumor.ab.expression.sort.cbind$wilcoxon.p.value,
  method = 'bonferroni'
  );

par(mfrow=c(1,3)); # Make three P-Values next to one another for comparison.

hist(tumor.ab.expression.sort.cbind$wilcoxon.p.value,
     main = 'Histogram of Wilcoxon\nRank Sum Test P-Values',
     freq = TRUE,
     xlab = 'P-Value',
     ylab ='Counts',
     ylim = c(0, 500),
     las = 1
     );

hist(tumor.ab.expression.sort.cbind$wilcoxon.p.value.FDR,
     main = 'Histogram of Wilcoxon\nRank Sum Test P-Values\nAdjusted by FDR',
     freq = TRUE,
     xlab = 'P-Value (FDR)',
     ylab = 'Counts',
     ylim = c(0, 500),
     las = 1
     );

hist(tumor.ab.expression.sort.cbind$wilcoxon.p.value.bonferroni,
     main = 'Histogram of Wilcoxon\nRank Sum Test P-Values\nAdjusted by Bonferroni',
     freq = TRUE,
     xlab = 'P-Value (Bonferroni)',
     ylab = 'Counts',
     ylim = c(0, 500),
     las = 1
     );

par(mfrow=c(1,1)); # Resets the mfrow parameter.


# Q5 Answers:
# Perform permutation-based test.

### Function calc.perm.median ######################################################################
# Description: Calculate median for randomly selected gene expression values from three samples.
# Input variable:
#   gene.expr: A numeric vector to represent gene expression values in all 12 samples.
# Output variable:
#   median(perm.selected): median gene expression.
calc.perm.median <- function(gene.expr){
  perm.selected <- sample(as.numeric(gene.expr),
                          3,
                          replace = FALSE
                          ); # Randomly selected three values from permutation.
  return(median(perm.selected)); # Calculate the median.
  }

# Calculate P-Values from permutation.
### Function run.permutation.pval.compare.exp.obs ##################################################
# Description: Calculate P-Values from permutation by comparing observed and expected values.
# Input variable:
#   gene.expr: A numeric vector to represent gene expression values in all 12 samples.
# Output variable:
#   perm.output: a list consisting of permutation P-Value, the observed gene expression value, and a vector of 1,000 expected expression values for each gene.
run.permutation.pval.compare.exp.obs <- function(gene.expr){
  tumor.a.gene <- as.numeric(gene.expr[tumor.a.patients]);
  tumor.b.gene <- as.numeric(gene.expr[tumor.b.patients]);
  
  tumor.a.obs <- median(tumor.a.gene); # Observed value: median gene expression in three Tumor subtype A samples.
  tumor.ab.obs <- median(c(tumor.a.gene, tumor.b.gene)); # Overall median gene expression among all 12 samples.
  
  # Generate 1,000 expected gene expression values, each of which is
  # a median from randomly selected three gene expression values out of all 12 samples.
  ### Function replicate #############################################################################
  # Description: Use replicate() function to repeatedly evaluate the output of calc.perm.median() function.
  # Input variable:
  #   n: the number of replications (1,000 times).
  #   expression: calc.perm.median(c(tumor.a.gene, tumor.b.gene)), which repeatedly evaluates overall median expression for each gene from random sampling.
  # Output variable:
  #   tumor.ab.exp: a numeric vector of 1,000 expected expression values for each gene.
  tumor.ab.exp <- replicate(n=1000,
                            calc.perm.median(c(tumor.a.gene, tumor.b.gene))
                            );
  
  # Perform two-sided permutation test: 
  # Calculate the number of times in which
  # abs(tumor.ab.exp - tumor.ab.obs), the difference between expected gene expression and overall median,
  # is greater than
  # abs(tumor.a.obs - tumor.ab.obs), the difference between observed gene expression and overall median.
  # Use abs() function to implement two-sided test.
  perm.p.value <- sum(abs(tumor.ab.exp - tumor.ab.obs) > abs(tumor.a.obs - tumor.ab.obs))/1000;
  perm.output <- list(perm.p.value = perm.p.value,
                      tumor.a.obs = tumor.a.obs,
                      tumor.ab.exp = tumor.ab.exp
                      );
  return(perm.output);
  }

# Run run.permutation.pval.compare.exp.obs() function for every row (gene).
run.permutation.pval.compare.exp.obs.out <- apply(tumor.ab.expression.sort.cbind,
                                                  1,
                                                  run.permutation.pval.compare.exp.obs # Invoke run.permutation.pval.compare.exp.obs() function.
                                                  );

# Unlist to generate matrix.
run.permutation.pval.compare.exp.obs.out.sapply <- sapply(run.permutation.pval.compare.exp.obs.out,
                                                          unlist
                                                          );

# Transpose and convert to data frame, in which row names are 500 gene IDs and column names include one permutation P-Value, one observed gene expression value, and 1,000 expected expression values for each gene.
run.permutation.pval.compare.exp.obs.out.sapply.df <- as.data.frame(t(run.permutation.pval.compare.exp.obs.out.sapply));

# Histogram of permutation P-Values.
hist(run.permutation.pval.compare.exp.obs.out.sapply.df$perm.p.value,
     main = 'Histogram of Permutation Test P-Values\n(Comparing Exp to Obs)',
     freq = TRUE,
     xlab = 'P-Value',
     ylab = 'Counts',
     ylim = c(0, 100),
     las = 1,
     breaks = seq(0, 1, by = 0.05) # Display X axis labels at an increment of 0.05 between 0 and 1.
     );

# Generate FDR-adjusted P-Values.
run.permutation.pval.compare.exp.obs.out.sapply.df$perm.p.value.fdr <- p.adjust(
  run.permutation.pval.compare.exp.obs.out.sapply.df$perm.p.value,
  method = 'fdr'
  );

# Write gene ID, observed median, expected median, p-value, and adjusted p-value to file in a tab-delimited format.
write.table(run.permutation.pval.compare.exp.obs.out.sapply.df,
            file = 'Shu.Tao.run.permutation.pval.compare.exp.obs.out.sapply.df.txt',
            quote = FALSE,
            sep = '\t',
            row.names = TRUE,
            col.names = NA
            );


# Q6 Answers:
# Use BoutrosLab.statistics.general for p-value extraction and BoutrosLab.plotting.general for plots.
# List all functions in each package.
ls("package:BoutrosLab.statistics.general");
ls("package:BoutrosLab.plotting.general");

### Function BoutrosLab.statistics.general::get.ttest.p ############################################
# Description: Use BoutrosLab.statistics.general::get.ttest.p() function
# in for-loop to calculate T Test P-Values.
# Input variable:
#   gene.expr: A numeric vector to represent gene expression values in all 12 samples.
#   tumor.a.index: Tumor subtype A indices.
#   tumor.b.index: Tumor subtype B indices.
# Output variable:
#   gene.ttest.pval: T Test P-Value for a gene.
# Write all T Test P-Values to a new column with name 'ttest.p.value.PCB'.
for (gene.counter in 1:nrow(tumor.ab.expression.sort.cbind)){
  gene.expr <- tumor.ab.expression.sort.cbind[gene.counter,];
  tumor.a.index <- c('Patient1', 'Patient2', 'Patient3');
  tumor.b.index <- c('Patient4', 'Patient5', 'Patient6',
                     'Patient7', 'Patient8', 'Patient9',
                     'Patient10', 'Patient11', 'Patient12'
                     );
  
  gene.ttest.pval <- BoutrosLab.statistics.general::get.ttest.p(gene.expr,
                                                                tumor.a.index,
                                                                tumor.b.index,
                                                                alternative = 'two.sided',
                                                                var.equal = FALSE,
                                                                paired = FALSE
                                                                );
  tumor.ab.expression.sort.cbind[gene.counter, "ttest.p.value.PCB"] <- gene.ttest.pval;
}

### Function BoutrosLab.plotting.general::create.histogram #########################################
# Description: Use BoutrosLab.plotting.general::create.histogram() function to create 
# a histogram of T Test P-Values.
# Input variable:
#   tumor.ab.expression.sort.cbind$ttest.p.value.PCB: A list of T Test P-Values.
# Output variable:
#   A histogram plot.
# Default type here is 'percent' while the 'hist' function used in R has the default type as 'counts'.
BoutrosLab.plotting.general::create.histogram(tumor.ab.expression.sort.cbind$ttest.p.value.PCB,
                                              main = 'Histogram of T Test P-Values (PCB)',
                                              xlab.label = 'P-Value',
                                              ylab.label = 'Percent',
                                              main.cex = 1,
                                              xlab.cex = 1,
                                              ylab.cex = 1,
                                              xaxis.cex = 1,
                                              yaxis.cex = 1,
                                              xaxis.tck = c(1,0),
                                              yaxis.tck = c(1,0)
                                              );

# Use 'apply' function to invoke calc.ttest.pval.PCB() function to calculate T Test P-Values.
### Function calc.ttest.pval.PCB ###################################################################
# Description: Use calc.ttest.pval.PCB() function to calculate T-Test P-Value for each gene.
# Input variable:
#   gene.expr: A numeric vector to represent gene expression in all 12 samples.
# Output variable:
#   gene.ttest.pval.PCB: T Test P-Value for a gene between two tumor subtypes.
calc.ttest.pval.PCB <- function(gene.expr) {
  tumor.a.index <- c('Patient1', 'Patient2', 'Patient3');
  tumor.b.index <- c('Patient4', 'Patient5', 'Patient6',
                     'Patient7', 'Patient8', 'Patient9',
                     'Patient10', 'Patient11', 'Patient12'
                     );
  
  gene.ttest.pval.PCB <- BoutrosLab.statistics.general::get.ttest.p(gene.expr,
                                                                    tumor.a.index,
                                                                    tumor.b.index,
                                                                    alternative = 'two.sided',
                                                                    var.equal = FALSE,
                                                                    paired = FALSE
                                                                    );
  return(gene.ttest.pval.PCB);
  }

# Use apply() function to invoke calc.ttest.pval() function to calculate T Test P-Value
# for each row (gene) of the input data frame tumor.ab.expression.sort.cbind.
ttest.p.value.by.wrapper.PCB <- apply(tumor.ab.expression.sort.cbind,
                                      1,
                                      calc.ttest.pval.PCB # Invoke calc.ttest.pval.PCB() function.
                                      );


### Function BoutrosLab.statistics.general::get.utest.p.and.medianfoldchange #######################
# Description: Use BoutrosLab.statistics.general::get.utest.p.and.medianfoldchange() function
# in for-loop to calculate Wilcoxon Rank Sum Test P-Values.
# Input variable:
#   gene.expr: A numeric vector to represent gene expression values in all 12 samples.
#   tumor.a.index: Tumor subtype A indices.
#   tumor.b.index: Tumor subtype B indices.
# Output variable:
#   gene.wilcoxtest.pval: a vector of two components, Wilcoxon Rank Sum Test P-Value and fold change.
# Write all Wilcoxon Rank Sum Test P-Values to a new column with name 'wilcox.p.value.PCB'.
# Write all fold change to a new column with name 'FoldChange.PCB'.
for (gene.counter in 1:nrow(tumor.ab.expression.sort.cbind)) {
  gene.expr <- tumor.ab.expression.sort.cbind[gene.counter,];
  tumor.a.index <- c('Patient1', 'Patient2', 'Patient3');
  tumor.b.index <- c('Patient4', 'Patient5', 'Patient6',
                     'Patient7', 'Patient8', 'Patient9',
                     'Patient10', 'Patient11', 'Patient12'
                     );
  
  gene.wilcoxtest.pval <- BoutrosLab.statistics.general::get.utest.p.and.medianfoldchange(gene.expr,
                                                                                          tumor.a.index,
                                                                                          tumor.b.index,
                                                                                          paired = FALSE,
                                                                                          alternative = "two.sided",
                                                                                          logged = FALSE
                                                                                          );
  
  tumor.ab.expression.sort.cbind[gene.counter, "wilcox.p.value.PCB"] <- gene.wilcoxtest.pval[1];
  tumor.ab.expression.sort.cbind[gene.counter, "FoldChange.PCB"] <- gene.wilcoxtest.pval[2];
  }

### Function BoutrosLab.plotting.general::create.histogram #########################################
# Description: Use BoutrosLab.plotting.general::create.histogram() function to create 
# a histogram of Wilcoxon Rank Sum Test P-Values.
# Input variable:
#   tumor.ab.expression.sort.cbind$wilcox.p.value.PCB: A list of Wilcoxon Rank Sum Test P-Values.
# Output variable:
#   A histogram plot.
BoutrosLab.plotting.general::create.histogram(tumor.ab.expression.sort.cbind$wilcox.p.value.PCB,
                                              main = 'Histogram of Wilcoxon Rank Sum Test P-Values (PCB)',
                                              xlab.label = 'P-Value',
                                              ylab.label = 'Percent',
                                              main.cex = 1,
                                              xlab.cex = 1,
                                              ylab.cex = 1,
                                              xaxis.cex = 1,
                                              yaxis.cex = 1,
                                              xaxis.tck = c(1,0),
                                              yaxis.tck = c(1,0)
                                              );

# Use 'apply' function to invoke calc.wilcoxontest.pval.and.fc.PCB() function to calculate T Test P-Values.
### Function calc.wilcoxontest.pval.and.fc.PCB #####################################################
# Description: Use calc.wilcoxontest.pval.and.fc.PCB() function to calculate Wilcoxon Rank Sum Test
# P-Value for each gene.
# Input variable:
#   gene.expr: A numeric vector to represent gene expression in all 12 samples.
# Output variable:
#   output.pval.fc.PCB: a list consisting of Wilcoxon Rank Sum Test P-Value and fold change for the gene.
calc.wilcoxontest.pval.and.fc.PCB <- function(gene.expr) {
  tumor.a.index <- c('Patient1', 'Patient2', 'Patient3');
  tumor.b.index <- c('Patient4', 'Patient5', 'Patient6',
                     'Patient7', 'Patient8', 'Patient9',
                     'Patient10', 'Patient11', 'Patient12'
                     );
  
  gene.wilcoxtest.pval.fc.PCB <- BoutrosLab.statistics.general::get.utest.p.and.medianfoldchange(gene.expr,
                                                                                                 tumor.a.index,
                                                                                                 tumor.b.index,
                                                                                                 paired = FALSE,
                                                                                                 alternative = "two.sided",
                                                                                                 logged = FALSE
                                                                                                 );
  
  output.pval.fc.PCB <- list(pval = gene.wilcoxtest.pval.fc.PCB[1],
                             fc = gene.wilcoxtest.pval.fc.PCB[2]
                             );
  return(output.pval.fc.PCB);
  }

# Use apply() function to invoke calc.wilcoxontest.pval.and.fc.PCB() function to
# calculate Wilcoxon Rank Sum Test P-Value for each row (gene) of
# the input data frame tumor.ab.expression.sort.cbind.
wilcoxontest.pval.and.fc.by.wrapper.PCB <- apply(tumor.ab.expression.sort.cbind,
                                                 1,
                                                 calc.wilcoxontest.pval.and.fc.PCB
                                                 );

# Reformat result into data frame.
wilcoxontest.pval.and.fc.by.wrapper.PCB.sapply <- sapply(wilcoxontest.pval.and.fc.by.wrapper.PCB,
                                                         unlist
                                                         );
wilcoxontest.pval.and.fc.by.wrapper.PCB.sapply.df <- as.data.frame(t(wilcoxontest.pval.and.fc.by.wrapper.PCB.sapply));

# Wilcoxon Rank Sum Test P-Value and fold change calculated by using Boutros R packages in apply wrapper.
wilcoxontest.pval.by.wrapper.PCB <- wilcoxontest.pval.and.fc.by.wrapper.PCB.sapply.df$pval;
wilcoxontest.fc.by.wrapper.PCB <- wilcoxontest.pval.and.fc.by.wrapper.PCB.sapply.df$fc;


# Q7 Answers:
# Repeat Q5 analysis using a fold-change criterion.
# Perform permutation test to compare expected gene expression to
# the expression in Tumor subtype A using a fold-change (1.02) criterion.
# The rationale and R codes are extremely similar to those in Q5, except for implementing the fold change criterion.

### Function calc.perm.median ######################################################################
# Description: Calculate median for randomly selected gene expression values from three samples.
# Input variable:
#   gene.expr: A numeric vector to represent gene expression values in all 12 samples.
# Output variable:
#   median(perm.selected): median gene expression.
calc.perm.median <- function(gene.expr){
  perm.selected <- sample(as.numeric(gene.expr),
                          3,
                          replace = FALSE
                          ); # Randomly selected three values from permutation.
  return(median(perm.selected)); # Calculate the median.
  }

# Calculate P-Values from permutation.
### Function run.permutation.pval.using.exp.to.obs.foldchange ######################################
# Description: Calculate P-Values from permutation using expected-to-observed fold change criterion.
# Input variable:
#   gene.expr: A numeric vector to represent gene expression values in all 12 samples.
# Output variable:
#   perm.output: a list consisting of permutation P-Value, the observed gene expression value, and a vector of 1,000 expected expression values for each gene.
run.permutation.pval.using.exp.to.obs.foldchange <- function(gene.expr){
  tumor.a.gene <- as.numeric(gene.expr[tumor.a.patients]);
  tumor.b.gene <- as.numeric(gene.expr[tumor.b.patients]);
  
  tumor.a.obs <- median(tumor.a.gene); # Observed value: median gene expression in three Tumor subtype A samples.
  tumor.ab.obs <- median(c(tumor.a.gene, tumor.b.gene)); # Overall median gene expression among all 12 samples.
  
  # Generate 1,000 expected gene expression values, each of which is
  # a median from randomly selected three gene expression values out of all 12 samples.
  ### Function replicate #############################################################################
  # Description: Use replicate() function to repeatedly evaluate the output of calc.perm.median() function.
  # Input variable:
  #   n: the number of replications (1,000 times).
  #   expression: calc.perm.median(c(tumor.a.gene, tumor.b.gene)), which repeatedly evaluates overall median expression for each gene from random sampling.
  # Output variable:
  #   tumor.ab.exp: a numeric vector of 1,000 expected expression values for each gene.
  tumor.ab.exp <- replicate(n=1000,
                            calc.perm.median(c(tumor.a.gene, tumor.b.gene))
                            );
  
  # Perform two-sided permutation test: 
  # Calculate the number of times in which tumor.ab.exp/tumor.a.obs >= 1.02 or <= 1/1.02.
  # Use abs() function to implement two-sided test.
  perm.p.value <- sum( abs(log2(tumor.ab.exp/tumor.a.obs)) >= log2(1.02) )/1000; # Fold change >= 1.02 or <= 1/1.02.
  perm.output <- list(perm.p.value = perm.p.value,
                      tumor.a.obs = tumor.a.obs,
                      tumor.ab.exp = tumor.ab.exp
                      );
  return(perm.output);
  }

# Run run.permutation.pval.using.exp.to.obs.foldchange() function for every row (gene).
run.permutation.pval.using.exp.to.obs.foldchange.out <- apply(tumor.ab.expression.sort.cbind,
                                                              1,
                                                              run.permutation.pval.using.exp.to.obs.foldchange # Invoke run.permutation.pval.using.exp.to.obs.foldchange() function.
                                                              );

# Unlist to generate matrix.
run.permutation.pval.using.exp.to.obs.foldchange.out.sapply <- sapply(run.permutation.pval.using.exp.to.obs.foldchange.out, unlist);

# Transpose and convert to data frame, in which row names are 500 gene IDs and column names include one permutation P-Value, one observed gene expression value, and 1,000 expected expression values for each gene.
run.permutation.pval.using.exp.to.obs.foldchange.out.sapply.df <- as.data.frame(t(run.permutation.pval.using.exp.to.obs.foldchange.out.sapply));

# Histogram of permutation P-Values.
hist(run.permutation.pval.using.exp.to.obs.foldchange.out.sapply.df$perm.p.value,
     main = 'Histogram of Permutation Test P-Values\n(Using Exp-to-Obs Fold Change)',
     freq = TRUE,
     xlab = 'P-Value',
     ylab = 'Counts',
     ylim = c(0, 100),
     las = 1,
     breaks = seq(0, 1, by = 0.05)
     );

# Generate FDR-adjusted P-Values.
run.permutation.pval.using.exp.to.obs.foldchange.out.sapply.df$perm.p.value.fdr <- p.adjust(run.permutation.pval.using.exp.to.obs.foldchange.out.sapply.df$perm.p.value, method = 'fdr');

# Write gene ID, observed median, expected median, p-value, and adjusted p-value to file in a tab-delimited format.
write.table(run.permutation.pval.using.exp.to.obs.foldchange.out.sapply.df,
            file = 'Shu.Tao.run.permutation.pval.using.exp.to.obs.foldchange.out.sapply.df.txt',
            quote = FALSE,
            sep = '\t',
            row.names = TRUE,
            col.names = NA
            );

### Part I Summary #################################################################################
# In this study, microarray was used to measure the mRNA level of 500 genes from two tumor subtypes (A and B). A total number of 33 differentially expressed genes (DEGs) were identified by two-sample T test (P-Value < 0.05). As an alternative statistical test, Wilcoxon rank-sum test detected 24 DEGs (P-Value < 0.05). For both tests, no DEGs were left after P-Value adjustment with either FDR or Bonferroni procedure (adjusted P-Value < 0.05). Nonetheless, the difference between the FDR-based method[1] and the Bonferroni method[2] lies in that the former adjusts P-Values in a way that allows for certain degree of false discoveries while the latter is so conservative (strict) that it pushes the P-Values to an even higher end (less likely of a significant P-Value). In the last step, a permutation-based method that generated gene expression value from 1,000 random selections and compared such value to the observed expression value derived from three subtype A tumors resulted in 60 and 30 DEGs with unadjusted and FDR-adjusted P-Value < 0.05 respectively.
# References
# [1] Benjamini Y and Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57:289–300.
# [2] Neyman J and Pearson ES (1928). On the use and interpretation of certain test criteria for purposes of statistical inference. Biometrika, 20A:175–240.
### Part I Summary #################################################################################


### PartII Plotting.Q2.txt #########################################################################
# 2a. Create scatterplot using the 'cars' dataset provided in the R Datasets package.
### Function BoutrosLab.plotting.general::create.scatterplot #######################################
# Description: Use BoutrosLab.plotting.general::create.scatterplot() function to generate scatterplot.
# Input variable:
#   formula: the formula used to extract the x & y components from the data frame.
#   data: the data frame to plot.
# Output variable:
#   filename: the file name for the plot output; if not specified, the plot will be shown interactively.
BoutrosLab.plotting.general::create.scatterplot(
  formula = dist ~ speed,
  data = cars,
  main = 'Distance to Stop Based on Car Speed',
  main.just = 'center',
  main.cex = 1,
  axes.lwd = 2, # Width of axes lines.
  
  # Set the distance of the plot to the top, bottom, left and right margin.
  top.padding = 4,
  bottom.padding = 4,
  left.padding = 4,
  right.padding = 4,
  
  pch = 19, # Plotting character: 19 stands for solid circle.
  cex = 0.75, # Plotting character size.
  
  xlab.label = 'Speed (Miles/Hour)',
  xlab.cex = 1,
  xaxis.cex = 1,
  ylab.label = 'Distance (Feet)',
  ylab.cex = 1,
  yaxis.cex = 1,
  
  xaxis.tck = c(2, 0), # x axis tick sizes for bottom and top respectively.
  yaxis.tck = c(2, 0), # y axis tick sizes for left and right respectively.
  
  xaxis.fontface = 'plain',
  yaxis.fontface = 'plain'
  );

# 2b. Create heatmap using the 'Loblolly' dataset provided in the R Datasets package.
# Reshape the original 'Loblolly' dataset into wide format.
### Function reshape ###############################################################################
# Description: Use reshape() function to reformat a data frame with repeated measurements.
# Input variable:
#   data: the data frame to be reformatted.
#   idvar and v.names: the rows and columns in the resulting data frame.
#   timevar: a variable that identifies multiple records from the same group.
#   direction: the format of the resulting data frame.
# Output variable:
#   loblolly.reshape.wide: the resulting data frame.
loblolly.reshape.wide <- reshape(
  Loblolly,
  v.names = "height", # Columns corresponding to Heights in different Ages.
  idvar = "Seed", # Rows identifying each Seed.
  timevar = "age", # Differentiates multiple Heights from the same Seed.
  direction = "wide" # Specifies reshaping direction.
  );

# Assign 1st column (Seed number) of loblolly.reshape.wide as its row names.
rownames(loblolly.reshape.wide) <- loblolly.reshape.wide[,1];

# Remove 1st column.
loblolly.reshape.wide <- loblolly.reshape.wide[,-1];

# Create heatmap.
### Function BoutrosLab.plotting.general::create.heatmap ###########################################
# Description: Use BoutrosLab.plotting.general::create.heatmap() function to generate heatmap.
# Input variable:
#   x: the input data frame to be plotted.
# Output variable:
#   filename: the file name for the plot output; if not specified, the plot will be shown interactively.
BoutrosLab.plotting.general::create.heatmap(
  x = loblolly.reshape.wide,
  colourkey.cex = 1, # Adjusts colourkey size.
  colourkey.labels.at = seq(0, 65, 5), # Specifies the tick-positions on the colourkey.
  xaxis.lab = NA, # Set to NA to take the row names (e.g., 301, 303) of loblolly.reshape.wide by default as the x axis labels.
  yaxis.lab = sapply(strsplit(colnames(loblolly.reshape.wide), '.', fixed = TRUE), unlist)[2,],
  xlab.label = 'Seed',
  ylab.label = 'Height',
  xaxis.cex = 1.5,
  yaxis.cex = 1.5,
  xaxis.fontface = 'plain',
  yaxis.fontface = 'plain',
  top.padding = 1,
  bottom.padding = 1,
  left.padding = 1,
  right.padding = 1,
  axis.xlab.padding = 1, # Padding between axis of plot and x label.
  clustering.method = 'complete',
  rows.distance.method = 'jaccard',
  cols.distance.method = 'jaccard',
  stratified.clusters.cols = list(c(1:7), c(8:14)), # Specifies the column locations of the columns to be combined into a strata.
  grid.row = TRUE, # Turns on interior grid-lines on the row level.
  grid.col = TRUE, # Turns on interior grid-lines on the column level.
  row.colour = 'white', # Specifies the color of interior grid-lines on the row level as white.
  col.colour = c(rep('white', 7), 'red', rep('white', 7)), # Specifies the color of interior grid-lines on the column level as white. And, the column in the middle to divide two clusters was colored as red.
  row.lwd = 1, # Specifies 1 as the width of the row level grid-lines.
  col.lwd = c(rep(1,7), 6, rep(1,7)), # Specifies 1 as the width of the column level grid-lines. And, the column in the middle to divide two clusters was assigned a width of 6.
  colour.scheme = default.colours(5, palette.type = 'spiral.morning'), # Specifies discrete colour schemes for tree heights.
  total.colours = 6 # Set total number of colours to be one more than the colors in the colour scheme to account for missing colour.
  );

# 2c. Create scatterplot using the 'ChickWeight' dataset provided in the R Datasets package.
### Function BoutrosLab.plotting.general::create.scatterplot #######################################
# Description: Use BoutrosLab.plotting.general::create.scatterplot() function to generate scatterplot.
# Input variable:
#   formula: the formula used to extract the x & y components from the data frame.
#   data: the data frame to plot.
#   groups: the grouping variable in the data frame.
# Output variable:
#   filename: the file name for the plot output; if not specified, the plot will be shown interactively.
BoutrosLab.plotting.general::create.scatterplot(
  formula = weight ~ Time,
  data = ChickWeight,
  groups = Diet, # Diet 1 through 4.
  col = c('black', 'green', 'blue', 'red'), # Assign colors of black, green, blue, red to Diet 1 through 4.
  xlab.label = 'Time',
  ylab.label = 'Weight',
  axis.key.padding = 2,
  xaxis.tck = c(1,0),
  yaxis.tck = c(1,0),
  xaxis.cex = 1,
  yaxis.cex = 1,
  xat = unique(ChickWeight$Time),
  xgrid.at = unique(ChickWeight$Time),
  cex = 0.5,
  pch = 19,
  type = c("p", "smooth"), # The data in the scatterplot are represented by points ("p"), and also a smooth ("smooth") line fitted to the points.
  add.grid = TRUE,
  lty = 1,
  
  # Create legend to the right.
  key = list(
    text = list(
      lab = c('Diet 1', 'Diet 2', 'Diet 3', 'Diet 4'),
      col = c('black', 'green', 'blue', 'red')
      ),
    lines = list(
      col = c('black', 'green', 'blue', 'red')
      ),
    space = 'right'
    ),
  );


### PartII Plotting.Q3.txt #########################################################################
# This R script was used to re-generate the Q3.SampleOutput plot.

# Read input data.
Q3.SeqControl.data <- read.table(
  file = 'Q3_SeqControl_data',
  header = TRUE,
  sep = '\t'
  );

# Reorder data frame by yes votes.
Q3.SeqControl.data.OrderByYesVotes <- Q3.SeqControl.data[
  order(Q3.SeqControl.data$yes.votes,
        decreasing = TRUE
        ),
  ];

# Force sample name to number, and then to matrix in order to make heatmap.
Sample.Name <- as.matrix(as.numeric(Q3.SeqControl.data.OrderByYesVotes$CPCG));

# Assign color to sample, and make sample heatmap.
sample.color <- list(
  rosybrown1 = 'CPCG0003P',
  rosybrown4 = 'CPCG0005P',
  red = 'CPCG0007P',
  darkred = 'CPCG0040P',
  darkorange = 'CPCG0047P',
  gold = 'CPCG0063P',
  darkolivegreen3 = 'CPCG0098P',
  darkgreen = 'CPCG0102P',
  aquamarine = 'CPCG0103P',
  cyan4 = 'CPCG0123P',
  dodgerblue = 'CPCG0183P',
  darkblue = 'CPCG0184P'
  );

Sample.CPCG.heatmap <- BoutrosLab.plotting.general::create.heatmap(
  x = Sample.Name,
  clustering.method = "none",
  scale.data = FALSE,
  colour.scheme = names(unlist(sample.color)),
  total.colours = 13,
  force.grid.col = TRUE,
  grid.col = TRUE,
  col.lwd = 3,
  axes.lwd = 4,
  print.colour.key = FALSE,
  yaxis.tck = 0,
  height = 1,
  same.as.matrix = TRUE
  );

# Assign color to sample prep procedure (FFPE or Frozen), and make sample prep heatmap.
prep.color <- list(
  white = 'CPCG0003P',
  white = 'CPCG0005P',
  white = 'CPCG0007P',
  white = 'CPCG0040P',
  white = 'CPCG0047P',
  white = 'CPCG0063P',
  white = 'CPCG0098P',
  darkslategrey = 'CPCG0102P',
  darkslategrey = 'CPCG0103P',
  white = 'CPCG0123P',
  white = 'CPCG0183P',
  white = 'CPCG0184P'
  );

Sample.prep.heatmap <- BoutrosLab.plotting.general::create.heatmap(
  x = Sample.Name,
  clustering.method = "none",
  scale.data = FALSE,
  colour.scheme = names(unlist(prep.color)),
  total.colours = 13,
  force.grid.col = TRUE,
  grid.col = TRUE,
  col.lwd = 3,
  axes.lwd = 4,
  print.colour.key = FALSE,
  yaxis.tck = 0,
  height = 1,
  same.as.matrix = TRUE
  );

# Make heatmap to show base quality greater than zero.
Sample.BasesQualGreaterThanZero <- as.matrix(Q3.SeqControl.data.OrderByYesVotes$X..Bases...0.quality);

Sample.BasesQualGreaterThanZero.heatmap <- BoutrosLab.plotting.general::create.heatmap(
  x = Sample.BasesQualGreaterThanZero,
  clustering.method = "none",
  scale.data = FALSE,
  colour.scheme = c('white', 'darkorange'),
  force.grid.col = TRUE,
  grid.col = TRUE,
  col.lwd = 3,
  axes.lwd = 4,
  print.colour.key = FALSE,
  yaxis.tck = 0,
  height = 1,
  same.as.matrix = TRUE
  );

# Make heatmap to show unique start points.
Sample.UniqStPts <- as.matrix(Q3.SeqControl.data.OrderByYesVotes$Unique.start.points);

Sample.UniqStPts.heatmap <- BoutrosLab.plotting.general::create.heatmap(
  x = Sample.UniqStPts,
  clustering.method = "none",
  scale.data = FALSE,
  colour.scheme = c('white', 'darkblue'),
  force.grid.col = TRUE,
  grid.col = TRUE,
  col.lwd = 3,
  axes.lwd = 4,
  print.colour.key = FALSE,
  yaxis.tck = 0,
  height = 1,
  same.as.matrix = TRUE
  );

### Function BoutrosLab.plotting.general::scientific.notation ######################################
# Description: Use BoutrosLab.plotting.general::scientific.notation() function to reformat maximum and minimum unique start points using scientific notation.
# Input variable:
#   x: the input number to be reformatted.
#   digits: the number of decimal places to keep.
#   type: the return format.
# Output variable:
#   max.unique.start.points or min.unique.start.points: a list consisting of base and exponent of the input number.
max.unique.start.points <- BoutrosLab.plotting.general::scientific.notation(
  x = max(Q3.SeqControl.data.OrderByYesVotes$Unique.start.points),
  digits = 2,
  type = 'list'
  );
max.unique.start.points.base <- max.unique.start.points$base;
max.unique.start.points.exponent <- max.unique.start.points$exponent;

min.unique.start.points <- BoutrosLab.plotting.general::scientific.notation(
  x = min(Q3.SeqControl.data.OrderByYesVotes$Unique.start.points),
  digits = 2,
  type = 'list'
  );
min.unique.start.points.base <- min.unique.start.points$base;
min.unique.start.points.exponent <- min.unique.start.points$exponent;

# Max and min unique start points label in legend.
### Function bquote ################################################################################
# Description: Use bquote() function to allow partial substitution of variable with its value in expressions.
# Input variable:
#   expr: expression consisting of both variable and constant.
# Output variable:
#   max.unique.start.points.lab or min.unique.start.points.lab: scientific notation label in the plot legend.
max.unique.start.points.lab = bquote(.(max.unique.start.points.base)~"x 10"^~.(max.unique.start.points.exponent));
min.unique.start.points.lab = bquote(.(min.unique.start.points.base)~"x 10"^~.(min.unique.start.points.exponent));

# Make heatmap to show average reads start.
Sample.AvgReadsSt <- as.matrix(Q3.SeqControl.data.OrderByYesVotes$Average.reads.start);

Sample.AvgReadsSt.heatmap <- BoutrosLab.plotting.general::create.heatmap(
  x = Sample.AvgReadsSt,
  clustering.method = "none",
  scale.data = FALSE,
  colour.scheme = c('white', 'deeppink'),
  force.grid.col = TRUE,
  grid.col = TRUE,
  col.lwd = 3,
  axes.lwd = 4,
  print.colour.key = FALSE,
  yaxis.tck = 0,
  height = 1,
  same.as.matrix = TRUE
  );

# Make all legends.
all.legends.Q3.SeqControl.data <- list(
  
  # Sample legend
  legend = list(
    colours = c(
      'rosybrown1',
      'rosybrown4',
      'red',
      'darkred',
      'darkorange',
      'gold',
      'darkolivegreen3',
      'darkgreen',
      'aquamarine',
      'cyan4',
      'dodgerblue',
      'darkblue'
      ),
    labels = c(
      'CPCG0003P',
      'CPCG0005P',
      'CPCG0007P',
      'CPCG0040P',
      'CPCG0047P',
      'CPCG0063P',
      'CPCG0098P',
      'CPCG0102P',
      'CPCG0103P',
      'CPCG0123P',
      'CPCG0183P',
      'CPCG0184P'
      ),
    title = expression(underline('Sample')),
    continuous = FALSE
    ),
  
  # Sample preparation legend
  legend = list(
    colours = c(
      'white',
      'darkslategrey'
      ),
    labels = c(
      'Frozen',
      'FFPE'
      ),
    title = expression(underline('Sample preparation')),
    continuous = FALSE
    ),
  
  # Bases quality greater than zero legend
  legend = list(
    colours = c(
      'darkorange',
      'white'
      ),
    labels = c(
      round(max(Q3.SeqControl.data.OrderByYesVotes$X..Bases...0.quality), 1),
      round(min(Q3.SeqControl.data.OrderByYesVotes$X..Bases...0.quality), 1)
      ),
    title = expression(underline('% Base > 0 quality')),
    continuous = TRUE,
    cex = 2,
    pos.y = 0.25
    ),
  
  # Unique start points legend
  legend = list(
    colours = c(
      'darkblue',
      'white'
      ),
    labels = sapply(
      c(max.unique.start.points.lab,
        min.unique.start.points.lab
        ),
      as.expression
      ),
    title = expression(underline('Unique start points')),
    continuous = TRUE,
    cex = 2,
    pos.y = 0.25
    ),
  
  # Average reads start legend
  legend = list(
    colours = c(
      'deeppink',
      'white'
      ),
    labels = c(
      round(max(Q3.SeqControl.data.OrderByYesVotes$Average.reads.start),3),
      round(min(Q3.SeqControl.data.OrderByYesVotes$Average.reads.start),3)
      ),
    title = expression(underline('Average reads/start')),
    continuous = TRUE,
    cex = 2,
    pos.y = 0.25
    )
  
  );

### Function BoutrosLab.plotting.general::legend.grob ##############################################
# Description: Use BoutrosLab.plotting.general::legend.grob() function to takes a list and
# generates a grob representing one or more legends.
# Input variable:
#   legends: a list defining one or more legends.
# Output variable:
#   legends.combined: a grob representing one or more legends.
legends.combined <- BoutrosLab.plotting.general::legend.grob(
  legends = all.legends.Q3.SeqControl.data,
  title.just = 'left',
  title.cex = 2,
  label.cex = 2,
  title.fontface = 'plain',
  between.row = 2
  );

# Boxplot legend
yes.votes.frac.boxplot.legend <- list(
  colours = c(
    'grey',
    'black'
    ),
  labels = c(
    '< 50x',
    as.expression(substitute(Observed>='50x', list(Observed = ''))) # Changes >= to its mathematical symbol.
    ),
  title = expression(underline('Observed'))
  );

yes.votes.frac.boxplot.legend.grob <- BoutrosLab.plotting.general::legend.grob(
  legends = list(legend = yes.votes.frac.boxplot.legend),
  title.just = 'centre',
  title.cex = 3,
  label.cex = 3,
  size = 8,
  x = 1.32,
  y = 0.45
  );

# Set bar colors based on 'outcome' column.
bar.colours.one <- replace(as.vector(Q3.SeqControl.data.OrderByYesVotes$outcome),
                           which(1 == Q3.SeqControl.data.OrderByYesVotes$outcome),
                           'black'
                           );
bar.colours.one.zero <- replace((bar.colours.one),
                       which(0 == bar.colours.one),
                       'grey'
                       );

### Function BoutrosLab.plotting.general::create.barplot ###########################################
# Description: Use BoutrosLab.plotting.general::create.barplot() function to generate yes votes fraction barplot.
# Input variable:
#   formula: the formula used to extract the x & y components from the data frame.
#   data: the data frame to plot.
# Output variable:
#   yes.votes.frac.boxplot: barplot object.
yes.votes.frac.boxplot <- BoutrosLab.plotting.general::create.barplot(
  formula = yes.votes ~ rownames(Q3.SeqControl.data.OrderByYesVotes),
  data = Q3.SeqControl.data.OrderByYesVotes,
  sample.order = 'decreasing',
  border.col = NULL,
  yat = seq(0, 1, by = 0.25),
  ylimits = c(0,1.025),
  xlab.label = NULL,
  ylab.label = 'Fraction of yes votes',
  ylab.cex = 4.2,
  xaxis.cex = 0,
  xaxis.tck = 0,
  yaxis.cex = 4.2,
  yaxis.tck = c(2,0),
  axes.lwd = 5,
  abline.h = 0.50,
  abline.lty = 2,
  abline.lwd = 3,
  abline.col = 'grey',
  col = bar.colours.one.zero,
  border.lwd = 0,
  legend = list(inside = list(fun = yes.votes.frac.boxplot.legend.grob))
  );

# Make multipanel plot.
pdf('Shu.Tao.PartII.Q3.SampleOutput.pdf',
    width = 24,
    height = 16
    );

### Function BoutrosLab.plotting.general::create.multipanelplot ####################################
# Description: Use BoutrosLab.plotting.general::create.multipanelplot() function to merges multiple plots.
# Input variable:
#   plot.objects: a list of plot objects.
#   legend: legend for the plot.
# Output variable:
#   filename: the file name for the plot output; if not specified, the plot will be shown interactively.
BoutrosLab.plotting.general::create.multipanelplot(
  plot.objects = list(
    yes.votes.frac.boxplot,
    Sample.CPCG.heatmap,
    Sample.prep.heatmap,
    Sample.BasesQualGreaterThanZero.heatmap,
    Sample.UniqStPts.heatmap,
    Sample.AvgReadsSt.heatmap
    ),
  plot.objects.heights = c(1, 0.1, 0.1, 0.1, 0.1, 0.1),
  plot.objects.widths = 0.8,
  xlab.axis.padding = 0.1,
  bottom.padding = 6,
  ylab.axis.padding = 6,
  legend = list(right = list(fun = legends.combined)),
  right.legend.padding = 0.5,
  right.padding = 1
  );

dev.off();


### PartII Plotting.Q4.txt #########################################################################
# This R script was used to re-generate the Q4.SampleOutput plot.

# Read input data.
Q4.HetStudy.data <- read.table(
  file = 'Q4_HetStudy_data.txt',
  header = TRUE,
  sep = '\t',
  row.names = 1,
  );

# Reorder 28 samples.
sample.order <- c(
  "CPCG0001F0",
  "CPCG0003F0",
  "CPCG0006F0",
  "CPCG0020F0",
  "CPCG0042F0",
  "CPCG0099F0",
  "CPCG0099P1",
  "CPCG0102F0",
  "CPCG0102P1",
  "CPCG0102P2",
  "CPCG0103P7",
  "CPCG0103F0",
  "CPCG0103P2",
  "CPCG0103P1",
  "CPCG0103P4",
  "CPCG0103P3",
  "CPCG0103P8",
  "CPCG0103P5",
  "CPCG0103P6",
  "CPCG0183F0",
  "CPCG0183P2",
  "CPCG0183P1",
  "CPCG0183P3",
  "CPCG0184P3",
  "CPCG0184P1",
  "CPCG0184F0",
  "CPCG0184P2",
  "CPCG0184P4"
);

# Also have to reorder Baca, Weischenfeldt and Berger so that they appear from bottom to top in the publication heatmap.
Q4.HetStudy.data.Reorder <- Q4.HetStudy.data[,c(sample.order,
                                                'Baca',
                                                'Weischenfeldt',
                                                'Berger'
                                                )
                                             ];

# Extract column index for samples.
sample.col = grep('^CPCG',
                  colnames(Q4.HetStudy.data.Reorder),
                  ignore.case = FALSE,
                  fixed = FALSE
                  );

# Calculate mutation fraction for every chromosomal region.
for(row.number in 1:nrow(Q4.HetStudy.data.Reorder)){
  chr.region.all = as.numeric(Q4.HetStudy.data.Reorder[row.number, sample.col]); # Extract mutation status in all samples.
  chr.region.withMut = length(chr.region.all[chr.region.all != 0]); # Calculate the number of samples having mutations.
  Q4.HetStudy.data.Reorder[row.number, 'Fraction'] = chr.region.withMut / length(chr.region.all); # Calculate mutation fraction.
}

# Make mutation fraction barplot.
fraction.barplot <- BoutrosLab.plotting.general::create.barplot(
  formula = Fraction ~ rownames(Q4.HetStudy.data.Reorder),
  data = Q4.HetStudy.data.Reorder,
  xlab.label = NULL,
  ylab.label = 'Fraction',
  xaxis.lab = rep('', nrow(Q4.HetStudy.data.Reorder)),
  xaxis.tck = 0,
  yaxis.tck = c(1,0),
  ylimits = c(0,0.5),
  yat = seq(0, 0.5, by = 0.25),
  );

# Specifies the chromosomes to be plotted.
chr.to.plot <- c('chr1', 'chr2',
                 'chr3', 'chr4',
                 'chr5', 'chr6',
                 'chr7', 'chr8',
                 'chr9', 'chr10',
                 'chr11', 'chr12',
                 'chr13', 'chr14',
                 'chr15', 'chr16',
                 'chr17', 'chr18',
                 'chr19', 'chr20',
                 'chr21', 'chr22',
                 'chrX', 'chrY'
                 )

# Initialize a list with the size (number of chromosomes) needed.
chr.region.count <- vector(mode = "list",
                           length = length(chr.to.plot) # 24 chromosomes.
                           );

# Use for-loop to calculate the number of regions in each chromosome for placing grid column lines.
for (chr.num in chr.to.plot){
  chr.region.count[[chr.num]] = as.numeric(
    table(sapply(strsplit(rownames(Q4.HetStudy.data.Reorder),
                          split = ':',
                          fixed = TRUE
                          ), unlist)[1,])[chr.num]
    );
  }

# Calculate column lines to define chromosome boundaries.
chr.boundaries <- c(unlist(chr.region.count)[1] + 0.5, # Grid column line between chr1 and chr2.
                    sum(unlist(chr.region.count)[1:2]) + 0.5,
                    sum(unlist(chr.region.count)[1:3]) + 0.5,
                    sum(unlist(chr.region.count)[1:4]) + 0.5,
                    sum(unlist(chr.region.count)[1:5]) + 0.5,
                    sum(unlist(chr.region.count)[1:6]) + 0.5,
                    sum(unlist(chr.region.count)[1:7]) + 0.5,
                    sum(unlist(chr.region.count)[1:8]) + 0.5,
                    sum(unlist(chr.region.count)[1:9]) + 0.5,
                    sum(unlist(chr.region.count)[1:10]) + 0.5,
                    sum(unlist(chr.region.count)[1:11]) + 0.5,
                    sum(unlist(chr.region.count)[1:12]) + 0.5,
                    sum(unlist(chr.region.count)[1:13]) + 0.5,
                    sum(unlist(chr.region.count)[1:14]) + 0.5,
                    sum(unlist(chr.region.count)[1:15]) + 0.5,
                    sum(unlist(chr.region.count)[1:16]) + 0.5,
                    sum(unlist(chr.region.count)[1:17]) + 0.5,
                    sum(unlist(chr.region.count)[1:18]) + 0.5,
                    sum(unlist(chr.region.count)[1:19]) + 0.5,
                    sum(unlist(chr.region.count)[1:20]) + 0.5,
                    sum(unlist(chr.region.count)[1:21]) + 0.5,
                    sum(unlist(chr.region.count)[1:22]) + 0.5, # Grid column line between chr22 and chrX.
                    sum(unlist(chr.region.count)[1:23]) + 0.5 # Grid column line between chrX and chrY.
                    );

# Make publication data heatmap.
publication.heatmap <- BoutrosLab.plotting.general::create.heatmap(
  Q4.HetStudy.data.Reorder[,c('Baca', 'Weischenfeldt', 'Berger')],
  clustering.method = 'none',
  scale.data = FALSE,
  print.colour.key = FALSE, # Do not print colour key for publication data heatmap.
  colour.scheme = c('white', 'red', 'green', 'blue'),
  total.colours = 5,
  yaxis.tck = 0,
  grid.row = TRUE,
  grid.col = TRUE,
  force.grid.col = TRUE,
  col.lines = chr.boundaries # Add column lines to define chromosome boundaries in heatmap.
  );

# make main heatmap.
main.heatmap <- BoutrosLab.plotting.general::create.heatmap(
  Q4.HetStudy.data.Reorder[,sample.col],
  clustering.method = 'none',
  scale.data = FALSE,
  print.colour.key = FALSE,
  colour.scheme = c('white', 'cornflowerblue','darkolivegreen4','darkred'),
  total.colours = 5,
  xaxis.tck = 0,
  yaxis.tck = c(1,0),
  xaxis.fontface = 'plain',
  yaxis.fontface = 'plain',
  xaxis.cex = 1,
  yaxis.cex = 1,
  xaxis.lab = c(seq(1, 22, by = 1), 'X', 'Y'),
  
  # Calculate coordinates to place xaxis tick labels (1 through 22, X, and Y).
  xat = c(round(0 + unlist(chr.region.count)[1]/2),
          round(unlist(chr.region.count)[1] + unlist(chr.region.count)[2]/2),
          round(sum(unlist(chr.region.count)[1:2]) + unlist(chr.region.count)[3]/2),
          round(sum(unlist(chr.region.count)[1:3]) + unlist(chr.region.count)[4]/2),
          round(sum(unlist(chr.region.count)[1:4]) + unlist(chr.region.count)[5]/2),
          round(sum(unlist(chr.region.count)[1:5]) + unlist(chr.region.count)[6]/2),
          round(sum(unlist(chr.region.count)[1:6]) + unlist(chr.region.count)[7]/2),
          round(sum(unlist(chr.region.count)[1:7]) + unlist(chr.region.count)[8]/2),
          round(sum(unlist(chr.region.count)[1:8]) + unlist(chr.region.count)[9]/2),
          round(sum(unlist(chr.region.count)[1:9]) + unlist(chr.region.count)[10]/2),
          round(sum(unlist(chr.region.count)[1:10]) + unlist(chr.region.count)[11]/2),
          round(sum(unlist(chr.region.count)[1:11]) + unlist(chr.region.count)[12]/2),
          round(sum(unlist(chr.region.count)[1:12]) + unlist(chr.region.count)[13]/2),
          round(sum(unlist(chr.region.count)[1:13]) + unlist(chr.region.count)[14]/2),
          round(sum(unlist(chr.region.count)[1:14]) + unlist(chr.region.count)[15]/2),
          round(sum(unlist(chr.region.count)[1:15]) + unlist(chr.region.count)[16]/2),
          round(sum(unlist(chr.region.count)[1:16]) + unlist(chr.region.count)[17]/2),
          round(sum(unlist(chr.region.count)[1:17]) + unlist(chr.region.count)[18]/2),
          round(sum(unlist(chr.region.count)[1:18]) + unlist(chr.region.count)[19]/2),
          round(sum(unlist(chr.region.count)[1:19]) + unlist(chr.region.count)[20]/2),
          round(sum(unlist(chr.region.count)[1:20]) + unlist(chr.region.count)[21]/2),
          round(sum(unlist(chr.region.count)[1:21]) + unlist(chr.region.count)[22]/2),
          round(sum(unlist(chr.region.count)[1:22]) + unlist(chr.region.count)[23]/2),
          round(sum(unlist(chr.region.count)[1:23]) + unlist(chr.region.count)[24]/2)
          ),
  
  xaxis.rot = 0, # Do not rotate x axis tick labels.
  yaxis.lab = unname(sapply(colnames(Q4.HetStudy.data.Reorder)[1:28], substring, 9, 10)),
  grid.row = TRUE,
  row.lines = c(seq(1.5, 5.5, by = 1), 7.5, 10.5, 19.5, 23.5),
  grid.col = TRUE,
  force.grid.col = TRUE,
  col.lines = chr.boundaries # Add column lines to define chromosome boundaries in heatmap.
  );

# Make legends.
all.legends.Q4.HetStudy.data <- list(
  
  # Patient ID legend
  legend = list(
    colours = c('blue',
                'purple',
                'green',
                'orange',
                'yellow',
                'black',
                'wheat4',
                'green4',
                'grey',
                'red4'
                ),
    labels = c('CPCG0001',
               'CPCG0003',
               'CPCG0006',
               'CPCG0020',
               'CPCG0042',
               'CPCG0099',
               'CPCG0102',
               'CPCG0103',
               'CPCG0183',
               'CPCG0184'
               ),
    title = expression(underline('Patient ID')),
    continuous = FALSE
    ),
  
  # Cohort legend
  legend = list(
    colours = c('pink', 'royalblue'),
    labels = c('Sx', 'Bx'),
    title = expression(underline('Cohort')),
    continuous = FALSE
    ),
  
  # Gleason score legend
  legend = list(
    colours = c('yellow', 'orange', 'red'),
    labels = c('3+4', '4+3','4+4'),
    title = expression(underline('Gleason score')),
    continuous = FALSE
    ),
  
  # Tissue type legend
  legend = list(
    colours = c(colours()[557], colours()[532]),
    labels = c('FFPE', 'Frozen'),
    title = expression(underline('Tissue type')),
    continuous = FALSE
    ),
  
  # Publication data legend
  legend = list(
    colours = c('red', 'blue', 'green'),
    labels = c('Baca', 'Weischenfeldt','Berger'),
    title = expression(underline('Publication')),
    continuous = FALSE
    )
  
  );

# Combine all legends.
legends.combined <- BoutrosLab.plotting.general::legend.grob(
  legends = all.legends.Q4.HetStudy.data,
  title.just = 'left',
  title.cex = 2,
  label.cex = 2,
  title.fontface = 'plain',
  size = 3,
  between.row = 2 # Add space between legend components.
  );

# Draw covariates on the right using heatmap.
# Use number 1 through 17 to represent discrete sample status below in order to make the heatmap:
#  1 to 10: represent sample name CPCG0001 to CPCG0184;
#  11 and 12: represent cohort Bx and Sx;
#  13, 14 and 15: represent gleason score 3+4, 4+3 and 4+4;
#  16 and 17: represent tissue type Frozen and FFPE.
covariates.on.the.right <- data.frame(
  patient.id = c(1, 2, 3, 4, 5, rep(6, 2), rep(7, 3), rep(8, 9), rep(9, 4), rep(10, 5)),
  cohort = c(rep(11, 5), rep(12, 23)),
  gleason.score = c(rep(13, 2), 14, 13, 14, 13, rep(14, 2), rep(15, 2), rep(13, 5), 14, 13, rep(14, 3), rep(13, 3), 15, rep(14, 3), 13),
  tissue.type = c(rep(16, 6), 17, 16, rep(17, 3), 16, rep(17, 7), 16, rep(17, 5), 16, rep(17, 2)),
  row.names = colnames(Q4.HetStudy.data.Reorder[sample.col])
  );

# Assign colors.
patient.id.colorscheme <- c('blue', 'purple', 'green', 'orange', 'yellow', 'black', 'wheat4', 'green4', 'grey', 'red4');
cohort.colorscheme <- c('royalblue', 'pink');
gleason.score.colorscheme <- c('yellow', 'orange', 'red');
tissue.type.colorscheme <- c(colours()[532], colours()[557]);

# Calculate cell position to add text '+'.
plus.sign.CPCG0001F0 <- which(colnames(t(covariates.on.the.right)) == 'CPCG0001F0', arr.ind = TRUE); # 1
plus.sign.CPCG0103P3 <- which(colnames(t(covariates.on.the.right)) == 'CPCG0103P3', arr.ind = TRUE); # 16
plus.sign.CPCG0103P8 <- which(colnames(t(covariates.on.the.right)) == 'CPCG0103P8', arr.ind = TRUE); # 17
plus.sign.CPCG0103P5 <- which(colnames(t(covariates.on.the.right)) == 'CPCG0103P5', arr.ind = TRUE); # 18
plus.sign.CPCG0103P6 <- which(colnames(t(covariates.on.the.right)) == 'CPCG0103P6', arr.ind = TRUE); # 19

right.covariates.heatmap <- BoutrosLab.plotting.general::create.heatmap(
  x = t(as.matrix(covariates.on.the.right)),
  clustering.method = 'none',
  scale.data = FALSE,
  print.colour.key = FALSE,
  grid.row = TRUE,
  grid.col = TRUE,
  
  # Add colors.
  colour.scheme = c(patient.id.colorscheme,
                    cohort.colorscheme,
                    gleason.score.colorscheme,
                    tissue.type.colorscheme
                    ),
  
  xaxis.tck = 0,
  yaxis.tck = 0,
  axes.lwd = 1,
  total.colours = 18,
  
  # Specify cells to add text '+'.
  row.pos = c(plus.sign.CPCG0001F0,
              plus.sign.CPCG0103P3,
              plus.sign.CPCG0103P8,
              plus.sign.CPCG0103P5,
              plus.sign.CPCG0103P6
              ),
  col.pos = 3,
  cell.text = '+',
  text.col = 'black'
  );

# Bottom covariate.
bottom.covariate <- data.frame(
  value = 1:4,
  row.names = c('None', 'CTX', 'ITX', 'INV')
  );

bottom.covariate.heatmap <- BoutrosLab.plotting.general::create.heatmap(
  x = t(as.matrix(bottom.covariate)),
  clustering.method = 'none',
  print.colour.key = FALSE,
  colour.scheme = c('white', 'cornflowerblue','darkolivegreen4','darkred'),
  total.colours = 5,
  xaxis.tck = 0,
  yaxis.tck = 0,
  xat = 1:4,
  xaxis.lab = rownames(bottom.covariate),
  xaxis.rot = 0,
  col.lwd = 1,
  grid.col = TRUE,
  axes.lwd = 1
  );

# Make multipanel plot.
pdf('Shu.Tao.PartII.Q4.HetStudy.data.Reorder.pdf',
    width = 20,
    height = 16
    );
BoutrosLab.plotting.general::create.multipanelplot(
  plot.objects = list(fraction.barplot,
                      publication.heatmap,
                      main.heatmap,
                      right.covariates.heatmap,
                      bottom.covariate.heatmap
                      ),
  layout.height = 4,
  layout.width = 2,
  layout.skip = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE), # Skip three locations to place five plots into eight locations.
  plot.objects.heights = c(3, 2, 15, 1.2),
  plot.objects.widths = c(16, 1.5),
  xlab.axis.padding = 0.1,
  bottom.padding = 6,
  ylab.axis.padding = 6,
  legend = list(left = list(fun = legends.combined)),
  left.legend.padding = 0.5,
  x.spacing = -3.5 # Adjust horizontal spacing between each plot.
  );

dev.off();






### 20200224 Updates ###############################################################################

### BPG Version Plot of the Histogram on Page 6 ####################################################
# Log10-based histogram of T-Test P-Values
BoutrosLab.plotting.general::create.histogram(log10(tumor.ab.expression.sort.cbind$ttest.p.value.PCB),
                                              main = 'Histogram of T Test P-Values (PCB)',
                                              xlab.label = 'Log10-based P-Value',
                                              ylab.label = 'Percent',
                                              main.cex = 1,
                                              xlab.cex = 1,
                                              ylab.cex = 1,
                                              xaxis.cex = 1,
                                              yaxis.cex = 1,
                                              xaxis.tck = c(1,0),
                                              yaxis.tck = c(1,0)
                                              );

### BPG Version Plot of the Histogram on Page 10 ###################################################
# Log10-based histogram of Wilcoxon Rank Sum Test P-Values
BoutrosLab.plotting.general::create.histogram(log10(tumor.ab.expression.sort.cbind$wilcox.p.value.PCB),
                                              main = 'Histogram of Wilcoxon Rank Sum Test P-Values (PCB)',
                                              xlab.label = 'Log10-based P-Value',
                                              ylab.label = 'Percent',
                                              main.cex = 1,
                                              xlab.cex = 1,
                                              ylab.cex = 1,
                                              xaxis.cex = 1,
                                              yaxis.cex = 1,
                                              xaxis.tck = c(1,0),
                                              yaxis.tck = c(1,0)
                                              );

### BPG Version Plot of the scatterplot on Page 11 #################################################
# Comparison of T-Test and Wilcoxon Rank Sum Test
# Create data frame to store test P-Values and test type; use number 1 and 2 for test type so that they can be jittered at next step.
t.and.wilcox.p.values.pcb <- data.frame(
  P.Value = c(tumor.ab.expression.sort.cbind$ttest.p.value,
              tumor.ab.expression.sort.cbind$wilcoxon.p.value
              ),
  Test = c(rep(1, # T Test
               length(tumor.ab.expression.sort.cbind$ttest.p.value)
               ),
           rep(2, # Wilcoxon Rank Sum Test
               length(tumor.ab.expression.sort.cbind$wilcoxon.p.value)
               )
           )
  );

### Function jitter ################################################################################
# Description: Use jitter() function to jitter points so that they are not stacked on one another on the vertical.
# Input variable:
#   x: a numeric vector, i.e., a vector of numbers indicating statistical test type.
#   factor: factor of the amount of noise to be added to a numeric vector.
# Output variable:
#   t.and.wilcox.p.values.pcb$Test.jitter: a vector of jittered numbers around each statistical type.
t.and.wilcox.p.values.pcb$Test.jitter <- jitter(t.and.wilcox.p.values.pcb$Test,
                                                factor = 2
                                                );

BoutrosLab.plotting.general::create.scatterplot(P.Value ~ Test.jitter,
                                                data = t.and.wilcox.p.values.pcb,
                                                groups = Test,
                                                main = 'Comparison of P-Values (PCB)',
                                                pch = 19, # Use solid circle
                                                cex = 0.5, # Solid circle size
                                                main.cex = 1,
                                                xlab.label = 'Test',
                                                xlab.cex = 1,
                                                ylab.cex = 1,
                                                xaxis.cex = 1,
                                                yaxis.cex = 1,
                                                xaxis.tck = c(1,0),
                                                yaxis.tck = c(1,0),
                                                xat = c(1,2),
                                                yat = seq(0, 1, by = 0.2),
                                                xaxis.lab = c('T', 'Wilcoxon Rank Sum'),
                                                ylimits = c(0,1)
                                                );

######## BPG version of three P-Values plots on Page 16 ############################################
pcb.hist.ttest.pval <- BoutrosLab.plotting.general::create.histogram(tumor.ab.expression.sort.cbind$ttest.p.value.PCB,
                                                                     xlab.label = 'P-Value',
                                                                     ylab.label = 'Percent',
                                                                     main.cex = 0.8,
                                                                     xlab.cex = 0.6,
                                                                     ylab.cex = 0.6,
                                                                     xaxis.cex = 0.6,
                                                                     yaxis.cex = 0.6,
                                                                     xaxis.tck = c(1,0),
                                                                     yaxis.tck = c(1,0),
                                                                     ylimits = c(0, 100),
                                                                     xat = seq(0, 1, by = 0.2)
                                                                     );

pcb.hist.ttest.pval.fdr <- BoutrosLab.plotting.general::create.histogram(p.adjust(tumor.ab.expression.sort.cbind$ttest.p.value.PCB,
                                                                                  method = 'fdr'
                                                                                  ),
                                                                         xlab.label = 'P-Value (FDR)',
                                                                         ylab.label = 'Percent',
                                                                         main.cex = 0.8,
                                                                         xlab.cex = 0.6,
                                                                         ylab.cex = 0.6,
                                                                         xaxis.cex = 0.6,
                                                                         yaxis.cex = 0.6,
                                                                         xaxis.tck = c(1,0),
                                                                         yaxis.tck = c(1,0),
                                                                         ylimits = c(0, 100),
                                                                         xat = seq(0, 1, by = 0.2)
                                                                         );

pcb.hist.ttest.pval.bonferroni <- BoutrosLab.plotting.general::create.histogram(p.adjust(tumor.ab.expression.sort.cbind$ttest.p.value.PCB,
                                                                                         method = 'bonferroni'
                                                                                         ),
                                                                                xlab.label = 'P-Value (Bonferroni)',
                                                                                ylab.label = 'Percent',
                                                                                main.cex = 0.8,
                                                                                xlab.cex = 0.6,
                                                                                ylab.cex = 0.6,
                                                                                xaxis.cex = 0.6,
                                                                                yaxis.cex = 0.6,
                                                                                xaxis.tck = c(1,0),
                                                                                yaxis.tck = c(1,0),
                                                                                ylimits = c(0, 100),
                                                                                xat = seq(0, 1, by = 0.2)
                                                                                );

BoutrosLab.plotting.general::create.multipanelplot(plot.objects = list(pcb.hist.ttest.pval,
                                                                       pcb.hist.ttest.pval.fdr,
                                                                       pcb.hist.ttest.pval.bonferroni
                                                                       ),
                                                   main = 'Histogram of T Test Raw, FDR-, and Bonferroni-Adjusted P-Values (PCB)',
                                                   main.cex = 0.8,
                                                   layout.height = 1, # Puts three histograms in one row in the final multipanelplot.
                                                   layout.width = 3, # Places three histograms in three positions of one row in the final multipanelplot.
                                                   x.spacing = 1,
                                                   ylab.axis.padding = 1,
                                                   width = 22 # Width of the final multipanelplot.
                                                   );

######## BPG version of three P-Values plots on Page 18 ############################################
pcb.hist.wilcox.pval <- BoutrosLab.plotting.general::create.histogram(tumor.ab.expression.sort.cbind$wilcox.p.value.PCB,
                                                                      xlab.label = 'P-Value',
                                                                      ylab.label = 'Percent',
                                                                      main.cex = 0.8,
                                                                      xlab.cex = 0.6,
                                                                      ylab.cex = 0.6,
                                                                      xaxis.cex = 0.6,
                                                                      yaxis.cex = 0.6,
                                                                      xaxis.tck = c(1,0),
                                                                      yaxis.tck = c(1,0),
                                                                      ylimits = c(0, 100),
                                                                      xat = seq(0, 1, by = 0.2)
                                                                      );

pcb.hist.wilcox.pval.fdr <- BoutrosLab.plotting.general::create.histogram(p.adjust(tumor.ab.expression.sort.cbind$wilcox.p.value.PCB,
                                                                                   method = 'fdr'
                                                                                   ),
                                                                          xlab.label = 'P-Value (FDR)',
                                                                          ylab.label = 'Percent',
                                                                          main.cex = 0.8,
                                                                          xlab.cex = 0.6,
                                                                          ylab.cex = 0.6,
                                                                          xaxis.cex = 0.6,
                                                                          yaxis.cex = 0.6,
                                                                          xaxis.tck = c(1,0),
                                                                          yaxis.tck = c(1,0),
                                                                          ylimits = c(0, 100),
                                                                          xat = seq(0, 1, by = 0.2)
                                                                          );

pcb.hist.wilcox.pval.bonferroni <- BoutrosLab.plotting.general::create.histogram(p.adjust(tumor.ab.expression.sort.cbind$wilcox.p.value.PCB,
                                                                                          method = 'bonferroni'
                                                                                          ),
                                                                                 xlab.label = 'P-Value (Bonferroni)',
                                                                                 ylab.label = 'Percent',
                                                                                 main.cex = 0.8,
                                                                                 xlab.cex = 0.6,
                                                                                 ylab.cex = 0.6,
                                                                                 xaxis.cex = 0.6,
                                                                                 yaxis.cex = 0.6,
                                                                                 xaxis.tck = c(1,0),
                                                                                 yaxis.tck = c(1,0),
                                                                                 ylimits = c(0, 100),
                                                                                 xat = seq(0, 1, by = 0.2)
                                                                                 );

BoutrosLab.plotting.general::create.multipanelplot(plot.objects = list(pcb.hist.wilcox.pval,
                                                                       pcb.hist.wilcox.pval.fdr,
                                                                       pcb.hist.wilcox.pval.bonferroni
                                                                       ),
                                                   main = 'Histogram of Wilcoxon Rank Sum Test\nRaw, FDR-, and Bonferroni-Adjusted P-Values (PCB)',
                                                   main.cex = 0.8,
                                                   layout.height = 1,
                                                   layout.width = 3,
                                                   x.spacing = 1,
                                                   ylab.axis.padding = 1,
                                                   width = 22
                                                   );

######## BPG version of three P-Values plots on Page 21 ############################################
BoutrosLab.plotting.general::create.histogram(run.permutation.pval.compare.exp.obs.out.sapply.df$perm.p.value,
                                              main = 'Histogram of Permutation Test P-Values\n(Comparing Exp to Obs)\n(PCB)',
                                              xlab.label = 'P-Value',
                                              ylab.label = 'Percent',
                                              main.cex = 1,
                                              xlab.cex = 1,
                                              ylab.cex = 1,
                                              xaxis.cex = 1,
                                              yaxis.cex = 1,
                                              xaxis.tck = c(1,0),
                                              yaxis.tck = c(1,0)
                                              );

######## BPG version of three P-Values plots on Page 34 ############################################
BoutrosLab.plotting.general::create.histogram(run.permutation.pval.using.exp.to.obs.foldchange.out.sapply.df$perm.p.value,
                                              main = 'Histogram of Permutation Test P-Values\n(Using Exp-to-Obs Fold Change)\n(PCB)',
                                              xlab.label = 'P-Value',
                                              ylab.label = 'Percent',
                                              main.cex = 1,
                                              xlab.cex = 1,
                                              ylab.cex = 1,
                                              xaxis.cex = 1,
                                              yaxis.cex = 1,
                                              xaxis.tck = c(1,0),
                                              yaxis.tck = c(1,0)
                                              );

#' 
#' 
## ----pressure, echo=TRUE------------------------------------------------------------------------------------------
### WRITE SESSION PROFILE TO FILE ##################################################################
BoutrosLab.utilities::save.session.profile(
  paste(date, 'BoutrosLab.Onboarding.Plot', 'Session.Info.txt', sep = '.')
  );

