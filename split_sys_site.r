############################################
######### input aa.txt #####################
##	seq1	ALVFGQM				############
##	seq2	ALVFGQM				############
##	seq3	ALVFGQM				############
#####     input gene.txt   #################
##	seq1	GCCCTTGTATTCGGCCAGATG	########
##	seq2	GCCCTAGTATTCGGCCAGATG	########
##	seq3	GCCCTTGTTTTCGGCCAGATG	########
############################################

### strsplit aa site
data1 <- read.table(file="aa.txt", header=F,row.name=1)
class(data1[,1])
aa <- as.character(data1[,1])
aaa <- strsplit(aa, "")  
aa_site  <- t(as.data.frame(aaa))
colnames(aa_site) <- c(1:ncol(aa_site))
rownames(aa_site) <- rownames(data1) 

### select synonymous position
table(aa_site[,1])
aa_num <- NULL
for(i in 1:ncol(aa_site)){
	aa_num[i] <- max(table(aa_site[,i]))
}
aa_num
select_aa_site <- which(aa_num==nrow(aa_site))
gene_site <- c((select_aa_site*3-2),(select_aa_site*3-1),select_aa_site*3)
select_site <- sort(as.vector(gene_site))

### select synonymous codon sites
data2 <- read.table(file="gene.txt", header=F,row.name=1)
class(data2[,1])
gene <- as.character(data2[,1])
gene1 <- strsplit(gene, "")  
gene2  <- t(as.data.frame(gene1))
colnames(gene2) <- c(1:ncol(gene2))
rownames(gene2) <- rownames(data2)
sys_site <- gene2[,select_site]
nonsys_site <- gene2[,-select_site]

### output synonymous codon sites
select_sys_site <- list()
for(i in 1:nrow(sys_site)){
	select_sys_site[i] <- paste(as.character(sys_site[i,1:ncol(sys_site)]),sep="",collapse="")
}
output_sys_site  <- t(as.data.frame(select_sys_site))
rownames(output_sys_site) <- rownames(sys_site)
write.table(output_sys_site, "sys_site.txt")

### output nonsynonymous codon sites
select_nonsys_site <- list()
for(i in 1:nrow(nonsys_site)){
	select_nonsys_site[i] <- paste(as.character(nonsys_site[i,1:ncol(nonsys_site)]),sep="",collapse="")
}
#select_nonsys_site
output_nonsys_site <- t(as.data.frame(select_nonsys_site))
rownames(output_nonsys_site) <- rownames(nonsys_site)
#output_nonsys_site
write.table(output_nonsys_site, "nonsys_site.txt")