library(plyr)
library(dada2)
library(seqinr)
library(digest)
ncpus<-as.numeric(commandArgs()[7])

# 7. Dereplicate reads
libs<-expand.grid(dir=c("FR","RF"),lib=c("R1","R2"),stringsAsFactors=F)
derep<-dlply(libs,.(lib),function(x) derepFastq(list.files("truncated",pattern=paste0(".*",unique(x$lib),".fastq"),full.names=T)) )

# 8. Learn error
err<-llply(derep,learnErrors,multithread=ncpus)

# 9. Apply dada2 error correction algorithm
dadas<-dlply(libs,.(dir,lib),function(x) dada(derep[[x$lib]][grep(x$dir,names(derep[[x$lib]]))],err[[x$lib]],pool=T,multithread=ncpus))

# 10. Merge pair-end reads
mergers<-dlply(libs,.(dir),function(x) {
  tmp<-apply(x,1,paste,collapse=".")
  mergePairs(dadas[[tmp[1]]],derep[["R1"]][grep(unique(x$dir),names(derep[["R1"]]))],
             dadas[[tmp[2]]],derep[["R2"]][grep(unique(x$dir),names(derep[["R2"]]))],
             minOverlap=10, maxMismatch=2, trimOverhang=T) %>%
    setNames(.,sub("\\..*","",names(.)))
})

# 11. Create count tables
seqtab_pairs<-llply(mergers,makeSequenceTable)
colnames(seqtab_pairs[["RF"]])<-laply(colnames(seqtab_pairs[["RF"]]),function(x) c2s(rev(comp(s2c(x),forceToLower=F))))

# 12. Merge count tables for sequences in both orientations
seqtab<-mergeSequenceTables(tables=seqtab_pairs,repeats="sum") %>%
  .[,nchar(colnames(.))>=300]

# 13. Remove bimeras
seqtab_clean<-removeBimeraDenovo(seqtab,method="pooled",multithread=ncpus)

# 14. Export ASVs to fasta file and count table to TAB-separated file
asv<-data.frame(asv=laply(colnames(seqtab_clean),digest::sha1),seq=colnames(seqtab_clean),stringsAsFactors=F)
colnames(seqtab_clean)<-asv$asv
write(apply(asv,1,function(x) paste0(">",paste(x,collapse="\n"))),file="V4_SSU_example.ASV.fasta",ncolumns=1)
write.table(data.frame(t(seqtab_clean)) %>% rownames_to_column("ASV"),"V4_SSU_example.ASV_table.tsv",sep="\t",col.names=T,row.names=F,quote=F)
