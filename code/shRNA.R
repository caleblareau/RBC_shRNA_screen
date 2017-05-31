#' ---
#' title: "Analyze shRNA screen with mixed linear models"
#' author: "Jacob C Ulirsch"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#' Load libraries
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
library(reshape2)
library(dplyr)
library(broom)
library(lme4)
library(nlme)
library(qvalue)
library(BuenColors)

#' Read in count data
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
shRNA <- read_delim("../data/countsForJacob.csv",delim=",")
shRNA <- as.data.frame(shRNA)
# Remove mouse genes
mouse.genes <- c("Psmb2", "Mtor", "Arcn1", "Rps7", "Nhp2l1", "Eif5b", "Eif2s2", "Ran", "Rpl5", "Psma1", "Rps11", "Rps4x", "Rps15a", "Rps9", "Snrpe", "Psmd1")
shRNA <- shRNA %>%
  filter(!(Gene %in% mouse.genes))
shRNA.old <- shRNA

#' Read in annotation data
annot <- read_delim("../data/shRNAref.csv",delim=",")
annot[annot$type=="Added",]$type <-"positive control"
annot[annot$type=="CP0009",]$type <-"positive control"
annot[annot$type=="CP0008",]$type <-"negative control"
annot[annot$type=="Test",]$type <-"GWAS"
annot.gene <- annot %>% 
  dplyr::select(type,gene,snp) %>% 
  unique()

#' Function to convert from counts to counts per million mapped reads
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
toCPM <- function(x) {
  log2(x/sum(x)*1000000+1)
}

#' Convert to CPM and set initial CPMs to 0
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
shRNA[,4:21] <- apply(shRNA[,4:21],2,toCPM)
shRNA[,4:9] <- shRNA[,4:9] - shRNA[,4]
shRNA[,10:15] <-shRNA[,10:15] - shRNA[,10]
shRNA[,16:21] <- shRNA[,16:21] - shRNA[,16]

#' Transform data from wide to long and create day and  donor variables
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
shRNA.melt <- melt(shRNA)
shRNA.melt$day <- as.numeric(gsub("r.*d","",shRNA.melt$variable))
shRNA.melt$donor <- gsub("d.*","",shRNA.melt$variable)
shRNA.old.melt <- melt(shRNA.old)
shRNA.old.melt$day <- as.numeric(gsub("r.*d","",shRNA.old.melt$variable))
shRNA.old.melt$donor <- gsub("d.*","",shRNA.old.melt$variable)

#' Pick only shRNAs that are well represented at the initial time point
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = FALSE
good.shRNAs <- shRNA.old.melt %>% 
  filter(day == 4) %>% 
  group_by(id) %>% 
  summarize(avg=mean(value),min=min(value)) %>% 
  filter(avg > 2^5, min > 2^4) %>% 
  .$id
shRNA.melt <- shRNA.melt %>%
  filter(id %in% good.shRNAs)

#' Pick only genes with enough shRNAs (4)
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
good.genes <- shRNA %>%
  filter(id %in% good.shRNAs) %>%
  group_by(Gene) %>%
  summarize(count=n()) %>%
  filter(count >= 4) %>%
  .$Gene
shRNA.melt <- shRNA.melt %>%
  filter(Gene %in% good.genes)

#' Run regression models allowing for each shRNA to have a random intercept and slope
#' Compare null (without fixed effect of time) to full model (including fixed effect of time)
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Have to use this slow method to make sure everything converges
# Can't use REML because comparing different fixed effects 
optimizer <- lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
shRNA.MM <- shRNA.melt %>%
  group_by(Gene) %>%
  do(beta = lmer(value ~ 1 + day + (0 + day | id), data=., control=optimizer, REML=FALSE)@beta[2],
     pvalue = anova(lmer(value ~ 1 + day + (0 + day | id), data=., control=optimizer, REML=FALSE),lmer(value ~ 1 + (0 + day | id), data=., REML=FALSE))[,8][2])
shRNA.MM$beta <- unlist(shRNA.MM$beta)
shRNA.MM$pvalue <- unlist(shRNA.MM$pvalue)
shRNA.MM$qvalue <- qvalue(shRNA.MM$pvalue)$qvalues

#' Annotate each gene
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
shRNA.MM.annot <- merge(shRNA.MM,annot.gene,by.x="Gene",by.y="gene")
shRNA.MM.annot$Gene <- factor(shRNA.MM.annot$Gene,ordered=T,
                              levels=unique(shRNA.MM.annot[order(shRNA.MM.annot$type,shRNA.MM.annot$pvalue),]$Gene))
shRNA.MM.nc <- shRNA.MM.annot %>% 
  filter(type == "negative control")
shRNA.MM.pc <- shRNA.MM.annot %>% 
  filter(type == "positive control")
shRNA.MM.GWAS <- shRNA.MM.annot %>% 
  filter(type == "GWAS")

#' Plot GWAS genes which are significant
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
sig.genes <- shRNA.MM.GWAS %>% 
  filter(pvalue < 0.05, beta < -0.15 | beta > 0.12) %>%
  .$Gene %>%
  as.vector()
pdf("../media/shRNA_sigGWASGenes.pdf")
for (i in sig.genes) {
  shRNA.melt.temp <- shRNA.melt %>% 
    filter(Gene == i)
  shRNA.MM.annot.temp <- shRNA.MM.annot %>%
    filter(Gene == i )
  beta.temp <- paste0("beta = ",signif(shRNA.MM.annot.temp$beta,3))
  pvalue.temp <- paste0("pvalue = ",signif(shRNA.MM.annot.temp$pvalue,3))
  p <- ggplot(shRNA.melt.temp,aes(y=value,x=day,group=interaction(id,donor),color=id)) + geom_line() +scale_fill_brewer(palette = "Spectral") + ggtitle(i) + theme_bw() + ylim(c(-10,4.6)) + labs(x="in vitro culture days", y="log2 CPM") +
    annotate("text", label = beta.temp, x = 4, y = 4, size = 5, hjust = 0) + annotate("text", label = pvalue.temp, x = 4, y = 3, size = 5, hjust = 0)
  print(p)
  }
dev.off()

#' Plot negative controls
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
nc.genes <- shRNA.MM.nc %>% 
  .$Gene %>%
  as.vector()
pdf("../media/shRNA_ncGenes.pdf")
for (i in nc.genes) {
  shRNA.melt.temp <- shRNA.melt %>% 
    filter(Gene == i)
  shRNA.MM.annot.temp <- shRNA.MM.annot %>%
    filter(Gene == i )
  beta.temp <- paste0("beta = ",signif(shRNA.MM.annot.temp$beta,3))
  pvalue.temp <- paste0("pvalue = ",signif(shRNA.MM.annot.temp$pvalue,3))
  p <- ggplot(shRNA.melt.temp,aes(y=value,x=day,group=interaction(id,donor),color=id)) + geom_line() +scale_fill_brewer(palette = "Spectral") + ggtitle(i) + theme_bw() + ylim(c(-10,4.6)) + labs(x="in vitro culture days", y="log2 CPM") +
    annotate("text", label = beta.temp, x = 4, y = 4, size = 5, hjust = 0) + annotate("text", label = pvalue.temp, x = 4, y = 3, size = 5, hjust = 0)
  print(p)
}
dev.off()

#' Plot positive controls
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pc.genes <- shRNA.MM.pc %>% 
  .$Gene %>%
  as.vector()
pdf("../media/shRNA_pcGenes.pdf")
for (i in pc.genes) {
  shRNA.melt.temp <- shRNA.melt %>% 
    filter(Gene == i)
  shRNA.MM.annot.temp <- shRNA.MM.annot %>%
    filter(Gene == i )
  beta.temp <- paste0("beta = ",signif(shRNA.MM.annot.temp$beta,3))
  pvalue.temp <- paste0("pvalue = ",signif(shRNA.MM.annot.temp$pvalue,3))
  p <- ggplot(shRNA.melt.temp,aes(y=value,x=day,group=interaction(id,donor),color=id)) + geom_line() +scale_fill_brewer(palette = "Spectral") + ggtitle(i) + theme_bw() + ylim(c(-10,4.6)) + labs(x="in vitro culture days", y="log2 CPM") +
    annotate("text", label = beta.temp, x = 4, y = 4, size = 5, hjust = 0) + annotate("text", label = pvalue.temp, x = 4, y = 3, size = 5, hjust = 0)
  print(p)
}
dev.off()

#' Plot number of hits per GWAS region
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
shRNA.MM.GWAS.snp <- shRNA.MM.GWAS %>% 
  filter(pvalue < 0.05, beta < -0.15 | beta > 0.12) %>%
  group_by(snp) %>%
  summarize(count=n()) %>% 
  arrange(count)
shRNA.MM.GWAS.snp$snp <- factor(shRNA.MM.GWAS.snp$snp,ordered=T,levels=shRNA.MM.GWAS.snp$snp)
pdf("../media/shRNA_perSNP.pdf")
ggplot(shRNA.MM.GWAS.snp,aes(y=count,x=snp,group=snp,fill=snp)) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position="none")
dev.off()

#' Get empirical FDR by sampling (it takes a while but is quite nice to have)
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
perms=100
GWASprobs <- shRNA.melt %>% 
  filter(Gene %in% unique(shRNA.MM.GWAS$Gene)) %>% 
  group_by(Gene) %>% 
  summarize(count = n_distinct(shRNA)) %>% 
  dplyr::select(count) %>% 
  table()
out <- data.frame(matrix(nrow = perms, ncol=2))
names(out) <- c("pvalue","beta")
shRNA.melt.nc <- shRNA.melt %>% 
  filter(day == 4, donor == "r1") %>% 
  filter(Gene %in% nc.genes) %>%
  .$id %>%
  unique()
for (i in seq(1,perms,1)) {
  perm1 <- shRNA.melt.nc %>% 
    sample(.,sample(as.numeric(names(GWASprobs)), size=1, replace=TRUE, prob=as.vector(GWASprobs/sum(GWASprobs))))
  perm1.out <- shRNA.melt %>%
    filter(id %in% perm1) %>%
    do(beta = lmer(value ~ 1 + day + (0 + day | id), data=., control=optimizer, REML=FALSE)@beta[2],
       pvalue = anova(lmer(value ~ 1 + day + (0 + day | id), data=., control=optimizer, REML=FALSE),lmer(value ~ 1 + (0 + day | id), data=., REML=FALSE))[,8][2])
  out[i,] <- c(perm1.out$pvalue,perm1.out$beta)
}
out.num <- out %>% filter(pvalue < 0.05, beta < -0.15 | beta > 0.12) %>% count() %>% .$n
# Empirical FDR
out.num / perms
