library(cluster)
library(Biobase)
library(qvalue)
NO_REUSE = F

# try to reuse earlier-loaded data if possible
if (file.exists("diffExpr.P1e-3_C2.matrix.RData") && ! NO_REUSE) {
    print('RESTORING DATA FROM EARLIER ANALYSIS')
    load("diffExpr.P1e-3_C2.matrix.RData")
} else {
    print('Reading matrix file.')
    primary_data = read.table("diffExpr.P1e-3_C2.matrix", header=T, com='', sep="\t", row.names=1, check.names=F)
    primary_data = as.matrix(primary_data)
}
source("/Volumes/air_clone/Applications/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/R/heatmap.3.R")
source("/Volumes/air_clone/Applications/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/R/misc_rnaseq_funcs.R")
source("/Volumes/air_clone/Applications/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/R/pairs3.R")
data = primary_data
sample_types = colnames(data)
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
sample_type_list = list()
for (i in 1:nsamples) {
    sample_type_list[[sample_types[i]]] = sample_types[i]
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
initial_matrix = data # store before doing various data transformations
data = log2(data+1)
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix
write.table(data, file="diffExpr.P1e-3_C2.matrix.log2.dat", quote=F, sep='	');
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
write.table(sample_cor, file="diffExpr.P1e-3_C2.matrix.log2.sample_cor.dat", quote=F, sep='	')
sample_dist = dist(t(data), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')
pdf("diffExpr.P1e-3_C2.matrix.log2.sample_cor_matrix.pdf")
sample_cor_for_plot = sample_cor
heatmap.3(sample_cor_for_plot, dendrogram='both', Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col = greenred(75), scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, symbreaks=F, margins=c(10,10), cexCol=1, cexRow=1, cex.main=0.75, main=paste("sample correlation matrix
", "diffExpr.P1e-3_C2.matrix.log2") )
dev.off()
gene_cor = NULL
gene_dist = dist(data, method='euclidean')
if (nrow(data) <= 1) { message('Too few genes to generate heatmap'); quit(status=0); }
hc_genes = hclust(gene_dist, method='complete')
myheatcol = colorpanel(75, 'purple','black','yellow')
data = t(scale(t(data), scale=F)) # center rows, mean substracted
write.table(data, file="diffExpr.P1e-3_C2.matrix.log2.centered.dat", quote=F, sep='	');
# redo the sample clustering according to the centered distance values.
sample_dist = dist(t(data), method='euclidean')
# redo the sample clustering
hc_samples = hclust(sample_dist, method='complete')
# redo the gene clustering according to centered distance values.
gene_dist = dist(data, method='euclidean')
# redo the gene clustering
hc_genes = hclust(gene_dist, method='complete')
heatmap_data = data
pdf("diffExpr.P1e-3_C2.matrix.log2.centered.genes_vs_samples_heatmap.pdf")
heatmap.3(heatmap_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, margins=c(10,10), cex.main=0.75, main=paste("samples vs. features
", "diffExpr.P1e-3_C2.matrix.log2.centered" ) )
dev.off()
save(list=ls(all=TRUE), file="diffExpr.P1e-3_C2.matrix.RData")
