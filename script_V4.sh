#bin/bash

# 0. Set environmental variables
NCPUS=32
FWD=TAReuk454FWD1
RVS=TAReukREV3

# 1. Split raw reads by orientation
cutadapt --no-indels --trimmed-only --action=none -O 15 -e 0.3 -g file:primers.txt -G file:primers.txt -o {name1}-{name2}.V4_SSU_example.lib1_R1.fastq -p {name1}-{name2}.V4_SSU_example.lib1_R2.fastq V4_SSU_example.lib1_R1.fastq.gz V4_SSU_example.lib1_R2.fastq.gz > log_cutadapt.primer_split.txt

# 2. For each orientation, demultiplex paired reads and remove barcode sequence
while read fwd rvs
do
	cutadapt --no-indels --trimmed-only -O 8 -e 0.25 -g file:${fwd}.barcodes.txt -G file:${rvs}.barcodes.txt -o {name1}-{name2}.${fwd}-${rvs}.V4_SSU_example.lib1_R1.fastq -p {name1}-{name2}.${fwd}-${rvs}.V4_SSU_example.lib1_R2.fastq ${fwd}-${rvs}.V4_SSU_example.lib1_R1.fastq ${fwd}-${rvs}.V4_SSU_example.lib1_R2.fastq > log_cutadapt.demultiplex.${fwd}_R1.${rvs}_R2.txt
done < <(paste <(sed -n '/>/s/>//p' primers.txt) <(sed -n '/>/s/>//p' primers.txt | tac))

# 3. Assign sample name for each pair of barcode and each orientation
while read sam fwd rvs
do
	mv ${fwd}-${rvs}.${FWD}-${RVS}.V4_SSU_example.lib1_R1.fastq ${sam}.FR.V4_SSU_example.lib1_R1.fastq
	mv ${fwd}-${rvs}.${FWD}-${RVS}.V4_SSU_example.lib1_R2.fastq ${sam}.FR.V4_SSU_example.lib1_R2.fastq
	mv ${rvs}-${fwd}.${RVS}-${FWD}.V4_SSU_example.lib1_R1.fastq ${sam}.RF.V4_SSU_example.lib1_R1.fastq
	mv ${rvs}-${fwd}.${RVS}-${FWD}.V4_SSU_example.lib1_R2.fastq ${sam}.RF.V4_SSU_example.lib1_R2.fastq
done < samples.txt

# 4. Remove primers at 5'-ends
mkdir trimmed
while read sam x y
do
	while read dir fwd rvs
	do
		cutadapt -j ${NCPUS} --no-indels --trimmed-only -O 15 -e 0.25 -g ^${fwd} -G ^${rvs} -o trimmed/${sam}.${dir}.V4_SSU_example.lib1_R1.trimmed.fastq -p trimmed/${sam}.${dir}.V4_SSU_example.lib1_R2.trimmed.fastq ${sam}.${dir}.V4_SSU_example.lib1_R1.fastq ${sam}.${dir}.V4_SSU_example.lib1_R2.fastq > log_cutadapt.trim_primer.${sam}.${dir}.txt
	done < <(paste <(echo -e "FR\nRF") <(sed -n '/>/!p' primers.txt) <(sed -n '/>/!p' primers.txt | tac))
done < samples.txt

# 5. Compute number of reads at different expected error and length thresholds
mkdir trim_stat
for i in trimmed/*
do
	out=$(basename $i)
	vsearch --quiet --fastq_eestats2 $i --length_cutoffs 210,*,5 --ee_cutoffs 0.5,1,1.5,2 --output - | sed '1,/^---/d;s/^  *//;s/  */\t/g;s/\([0-9]*\)(\t*\([0-9\.]*\)%)/\1 \2/g' > trim_stat/${out%.*}.stat
done

# 6. Truncate reads at the selected length and expected errors (e.g. lowest expected-error then highest sequencing length, above 210 nt one strand and above 460 nt both strands, at which all samples keep at least 80% of their sequences)
mkdir truncated
echo -e "FR 1.5 250 1.5 230\nRF 1 265 2 210" > truncate.parameters.txt
while read sam x y
do
	while read dir fee fl ree rl
	do
		vsearch --quiet --fastq_filter trimmed/${sam}.${dir}.V4_SSU_example.lib1_R1.trimmed.fastq --fastq_trunclen ${fl} --fastq_maxee ${fee} --fastq_maxns 0 --fastqout tmp_trunc.${sam}.${dir}.V4_SSU_example.lib1_R1.fastq
		vsearch --quiet --fastq_filter trimmed/${sam}.${dir}.V4_SSU_example.lib1_R2.trimmed.fastq --fastq_trunclen ${rl} --fastq_maxee ${ree} --fastq_maxns 0 --fastqout tmp_trunc.${sam}.${dir}.V4_SSU_example.lib1_R2.fastq
		join <(sed -n '1~4{s/[\t ].*$//;p}' tmp_trunc.${sam}.${dir}.V4_SSU_example.lib1_R1.fastq | sort) <(sed -n '1~4{s/[\t ].*$//;p}' tmp_trunc.${sam}.${dir}.V4_SSU_example.lib1_R2.fastq | sort) | sed 's/^@//' > ${sam}.${dir}.trunc.accnos
		vsearch --quiet --fastx_getseqs tmp_trunc.${sam}.${dir}.V4_SSU_example.lib1_R1.fastq --labels ${sam}.${dir}.trunc.accnos --fastqout truncated/${sam}.${dir}.V4_SSU_example.lib1_R1.fastq
		vsearch --quiet --fastx_getseqs tmp_trunc.${sam}.${dir}.V4_SSU_example.lib1_R2.fastq --labels ${sam}.${dir}.trunc.accnos --fastqout truncated/${sam}.${dir}.V4_SSU_example.lib1_R2.fastq
	done < truncate.parameters.txt
done < samples.txt
rm *trunc.accnos tmp_trunc.*

# 7. to 14. in R (or open a R terminal and execute commands from script_dada2.R manualy)
Rscript --vanilla script_dada2.R $NCPUS

# 15. Assign ASVs to PR2 reference taxonomy
wget https://github.com/vaulot/pr2_database/releases/download/v4.13.0/pr2_version_4.13.0_16S_UTAX.fasta.gz
wget https://github.com/vaulot/pr2_database/releases/download/v4.13.0/pr2_version_4.13.0_18S_UTAX.fasta.gz
gunzip -ck pr2_version_4.13.0_18S_UTAX.fasta.gz pr2_version_4.13.0_16S_UTAX.fasta.gz | vsearch --makeudb_usearch - --output pr2_version_4.13.0_UTAX.udb
vsearch --no_progress --usearch_global V4_SSU_example.ASV.fasta --threads ${NCPUS} --db pr2_version_4.13.0_UTAX.udb --dbmask none --qmask none --rowlen 0 --notrunclabels --userfields query+id2+target --maxaccepts 0 --maxrejects 32 --top_hits_only --output_no_hits --id 0.6 --iddef 2 --userout V4_SSU_example.ASV.hits

# 16. Compute least-common ancestor for ASV with multiple best matches at consensus threshold of 60 %
CONS=60
lca() {
	parallel --recstart ">" --remove-rec-sep --pipe -k -N1 awk -v cons=$CONS -f lca_vsearch.awk
}
export -f lca
export CONS
sed 's/;\t/\t/;s/;$//;s/;tax=/\t/' V4_SSU_example.ASV.hits | sort -k 1,1 -k 4,4 | awk '$1 != p{printf ">"}{p=$1}1' | parallel -j ${NCPUS} --recstart ">" --pipe -k lca > V4_SSU_example.ASV.taxonomy

# 17. Convert ASV count table to biom format, add all taxonomy information as observation metadata
biom convert -i V4_SSU_example.ASV_table.tsv -o tmp.biom  --table-type="OTU table" --to-json
sed 's/NA;/NA(0);/g;s/)*; / /;s/);/|/g;s/(/##/g' V4_SSU_example.ASV.taxonomy | awk '{l=split($3,a,"|");tax="";boot="";for(i=1;i<=l;i++){split(a[i],b,"##");tax=tax""b[1]";";boot=boot""b[2]";"};print $1,$2,tax,boot,$4}' | sed 's/; / /g;s/ /\t/g' | cat <(echo -e "#OTUID\tsimilarity\ttaxonomy\tbootstrap\treference") - > V4_SSU_example.ASV.taxo
biom add-metadata -i tmp.biom -o V4_SSU_example.ASV_table.biom --output-as-json --observation-metadata-fp V4_SSU_example.ASV.taxo --sc-separated taxonomy,bootstrap
