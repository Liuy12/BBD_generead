library('BSgenome.Hsapiens.UCSC.hg19')
library(MutationalPatterns)
library(ggplot2)
library(NMF)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
sample_files <- read.delim(sample_files, header = F, stringsAsFactors = F)
sample_names <- sample_files$V1
vcf_files <- sample_files$V2

message(paste0('Loading vcf files as granges object', '/n'))
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
## count mutation type occurrence
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
## 96 mutational profiles
message(paste0('saving mutation signature figure', '/n'))
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
gp <-  plot_96_profile(mut_mat[,1:n], condensed = TRUE)
ggsave(paste0(output,'/mutation_signature.png'), gp)

## de novo mutation signature extracting using NMF
message(paste0('saving de novo signatures', '\n'))
nmf_res <- extract_signatures(mut_mat, rank = nsig, nrun = 50)
colnames(nmf_res$signatures) <- paste0('Signature ', 1:nsig)
rownames(nmf_res$contribution) <- paste0('Signature ', 1:nsig)
## plot extracted signatures
gp <- plot_96_profile(nmf_res$signatures, condensed = TRUE)
ggsave(paste0(output,'/denovo_mutation_signatuers.png'),gp)
## plot contribution heatmap
png('heatmap_contmat_classical.png',width = 800,height=480)
aheatmap(t(contmat_classic),color = colorRampPalette(c('blue','white','red'))(1000),scale = 'r1',labRow = NA)
dev.off()
write.table(nmf_res$contribution, paste0(output, '/mutsig_classic.txt'),row.names=T, col.names=T,quote=F, sep='\t')
