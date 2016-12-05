library(Matrix)
library(Rcpp)
library(RcppEigen)

Rcpp::sourceCpp("scott_sgdlogit.cpp")

raw_data = read.csv("master.csv", header = TRUE)

# Read specific columns from the data

kegg = as.character(raw_data[ , 1])
gene = as.character(raw_data[ , 2])

data = raw_data[ , 3:ncol(raw_data)]

N = ncol(data)
P = nrow(data)
M = rep(1, N) # Used in logistic regression

# For C++ code, it's easiest to use a sparse matrix
# TODO: might want to consider scaling the counts, either in the way DESeq2 does
# it or something similar.
X = Matrix(as.matrix(unname(data)), sparse = TRUE)

# Convert the column names to 0 = control, 1 = treatment
Ynames = names(data)
Y = stringi::stri_endswith(Ynames, fixed = "T") * 1

# Initial guess for the betas.
beta_init = rep(0, P)
# For the weighting in the skip scenario for a feature
# TODO: algorithm is very sensitive to this value
eta = 0.0002
# passes through the data
nIter = 100
# L1 regularization
# TODO: need cross-validation or other way to choose lambda.
# Right now, using manual value by eye (BAD practice)
lambda = 1200.0
# exponentially-weighted moving average factor
# TODO: could explore other values.
discount = 0.01


# X[i, j] == data[i, j] is non-negative, so this trims out the genes with no counts
X_nozero = X[which(rowSums(data) > 0), ]
# Run the algorithm
result = sparsesgd_logit(X_nozero, Y, M, eta, nIter, beta_init, lambda, discount)

cat(paste("How many genes are in the final model? ", sum(result$beta != 0)))
cat(paste("Which genes?", paste(gene[which(result$beta != 0)], collapse = "\n"), sep = "\n"))
plot(result$nll_tracker)







library(DESeq2)
# Get the counts
count_data_raw = raw_data[ , 3:ncol(raw_data)]

########################################
#  Format the data in the way DESeq2 expects
########################################
new_colnames = rep(NULL, N)
condition = rep(NULL, N)
control = 0
treatment = 0
for (i in 1:N) {
  base = c("control", "treatment")[Y[i] + 1]
  val = 0
  if (base == "control") {
    control = control + 1
    val = control
  } else {
    treatment = treatment + 1
    val = treatment
  }
  new_colnames[i] = paste0(base, val)
  condition[i] = base
}
count_data = count_data_raw
colnames(count_data) = new_colnames
rownames(count_data) = gene

col_data = data.frame(condition)
rownames(col_data) = new_colnames

# Run DESeq2
dds = DESeqDataSetFromMatrix(countData = count_data,
                             colData = col_data,
                             design = ~ condition)
dds = DESeq(dds)
res = results(dds)

# plot(sort(res$pvalue))
# points(sort(res$padj), col = "red")
# P - sum(is.na(res$pvalue))
# P - sum(is.na(res$padj))

cat(paste("How many genes in default DESeq2?", paste(sort(gene[which(res$padj < 0.1)]), collapse = "\n"), sep="\n")) # 50 genes here

# Not exactly sure what this plot does either.
# plotMA(res, main="DESeq2", ylim=c(-3,3))

# Benjamini-Hochberg by hand (DESeq2 does something similar...)
alpha = 0.1
sorted_pval = sort(res$pvalue, na.last = TRUE)
numNA = sum(is.na(res$pvalue))
bh = sapply(1:P, function(i){ sorted_pval[i] <= alpha * (i + 1) / (P - numNA)})
threshold = sorted_pval[max(which(bh))]

cat(paste("How many genes in customized Benjamini-Hochberg procedure",
      paste(sort(gene[which(res$pvalue < threshold)]), collapse="\n"), sep="\n")) # 4 genes here

# These two are equal -- throws out genes that were never observed section 1.5.3 of DESeq2 paper
print(paste(sum(rowSums(count_data) == 0),
sum(is.na(res$pvalue))))

# 111675 genes are thrown out by the "independent filtering" for having a low
# mean normalized count
sum(is.na(res$padj)) - sum(is.na(res$pvalue))

# Not entirely sure what this plot means.
# plot(metadata(res)$filterNumRej,
#      type = "b", ylab = "number of rejections",
#      xlab = "quantiles of filter")
# lines(metadata(res)$lo.fit, col = "red")
# abline(v = metadata(res)$filterTheta)

#nofilter
# If you turn off the filtering, only two of the adjusted p-values are less than 0.1
# where the adjustment basically just accounts for Benjamini-Hochberg, I think
resNoFilt <- results(dds, independentFiltering = FALSE)
addmargins(table(filtering = (res$padj < .1),
noFiltering = (resNoFilt$padj < .1)))

cat(paste("How many genes in customized Benjamini-Hochberg procedure",
paste(sort(gene[which(resNoFilt$padj < 0.1)]), collapse="\n"),sep="\n")) # 2 genes here
