seqs <- readDNAStringSet("combined.fasta")
alignment <- msa(seqs)
print(alignment)

library(Biostrings)
library(seqLogo)

pfm <- consensusMatrix(DNAStringSet(seqs))[1:4,]
pfm <- prop.table(pfm,2)

seqLogo(pwm)

library(ape)
dna <- read.dna("aligned.fasta", format = "fasta")
dna

dist_matrix <- dist.dna(dna)
tree <- nj(dist_matrix)
plot(tree)

seqs<- readDNAStringSet("combined.fasta")
alignment <- msa(seqs)
aligned <- msaConvert(alignment,type = "ape::DNAbin")
dist_matrix <- dist.dna(aligned)
tree <- nj(dist_matrix)
plot(tree)

