library(Matrix)
library(Rcpp)
library(RcppEigen)

Rcpp::sourceCpp("scott_sgdlogit.cpp")

raw_data = read.csv("master.csv", header = TRUE)

kegg = as.character(raw_data[ , 1])
gene = as.character(raw_data[ , 2])

data = raw_data[ , 3:ncol(raw_data)]

N = ncol(data)
P = nrow(data)
M = rep(1, N)


X = Matrix(as.matrix(unname(data)), sparse = TRUE)
Ynames = names(data)
Y = stringi::stri_endswith(Ynames, fixed = "T") * 1

beta_init = rep(0, P)
# For the weighting in the skip scenario for a feature
# TODO: algorithm is very sensitive to this value
eta = 0.0002
# passes through the data
nIter = 100
# L1 regularization
# TODO: need CV or other way to choose lambda
lambda = 1200.0
# exponentially-weighted moving average
# TODO: could explore other values.
discount = 0.01

# X[i, j] == data[i, j] is non-negative
X_nozero = X[which(rowSums(data) > 0), ]
result = sparsesgd_logit(X_nozero, Y, M, eta, nIter, beta_init, lambda, discount)
#(MapMatd X, VectorXd Y, VectorXd M, double eta, int npass, VectorXd beta0, double lambda=1.0, double discount = 0.01)

sum(result$beta != 0)
gene[which(result$beta != 0)]
plot(result$nll_tracker)




library(DESeq2)
library(BiocParallel)
register(MulticoreParam(4))
count_data_raw = raw_data[ , 3:ncol(raw_data)]

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


dds = DESeqDataSetFromMatrix(countData = count_data,
                             colData = col_data,
                             design = ~ condition)
dds
dds = DESeq(dds, parallel = TRUE)
res = results(dds, parallel = TRUE)

plot(sort(res$pvalue))
points(sort(res$padj), col = "red")
P - sum(is.na(res$pvalue))
P - sum(is.na(res$padj))



plotMA(res, main="DESeq2", ylim=c(-3,3))

sort(gene[which(res$padj < 0.1)]) # 50 genes here

# Benjamini-Hochberg instead
alpha = 0.1
sorted_pval = sort(res$pvalue, na.last = TRUE)
numNA = sum(is.na(res$pvalue))
bh = sapply(1:P, function(i){ sorted_pval[i] <= alpha * (i + 1) / (P - numNA)})
threshold = sorted_pval[max(which(bh))]
sort(gene[which(res$pvalue < threshold)]) # 4 genes here

# These two are equal -- throws out genes that were never observed \S1.5.3 of paper
sum(rowSums(count_data) == 0)
sum(is.na(res$pvalue))

# 111675 genes are thrown out by the "independent filtering" for having a low
# mean normalized count
sum(is.na(res$padj)) - sum(is.na(res$pvalue))

plot(metadata(res)$filterNumRej,
     type = "b", ylab = "number of rejections",
     xlab = "quantiles of filter")
lines(metadata(res)$lo.fit, col = "red")
abline(v = metadata(res)$filterTheta)

#nofilter
# If you turn off the filtering, only two of the adjusted p-values are less than 0.1
# where the adjustment basically just accounts for Benjamini-Hochberg, I think
resNoFilt <- results(dds, independentFiltering = FALSE)
addmargins(table(filtering = (res$padj < .1),
                 noFiltering = (resNoFilt$padj < .1)))
sort(gene[which(resNoFilt$padj < 0.1)]) # only 2 genes here




