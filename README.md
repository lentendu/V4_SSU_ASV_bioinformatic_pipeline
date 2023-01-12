# Bioinformatic processing for eukaryotic 18S V4 Illumina 2x300nt raw sequencing data
 **Supplementary online data attached to the book chapter:**

*Lentendu, G., Lara, E., Geisen, S.*, 2023. **[Metabarcoding Approaches For Soil Eukaryotes, Protists and Microfauna](https://doi.org/10.1007/978-1-0716-2871-3_1)**, in: Martin, F., Uroz, S. (Eds.), Microbial Environmental Genomics (MEG), Methods in Molecular Biology, pp. 1–16. Springer, New York, NY

This repository provide the bash and R scripts as well as example data to produce a community matrix with taxonomic assignment out of raw Illumina reads, following all steps described in book chapter section 3.3.1. for case b.

The scripts were tested on a Debian GNU/Linux v.9.7 server with 32 CPUs and 124G of RAM.


## Dependencies:

 - in bash
   * cutadapt (https://cutadapt.readthedocs.io)
   * vsearch (https://github.com/torognes/vsearch)
   * GNU parallel (https://www.gnu.org/software/parallel/)
   * biom-format (https://biom-format.org/index.html)

 - in R
   * R (https://cran.r-project.org/)
   * plyr (https://cran.r-project.org/web/packages/plyr/index.html)
   * dada2 (https://benjjneb.github.io/dada2/index.html)
   * seqinr (https://cran.r-project.org/web/packages/seqinr/index.html)
   * digest (https://cran.r-project.org/web/packages/digest/index.html)


## Usage instruction:

Clone the repository:
```
git clone https://github.com/lentendu/V4_SSU_ASV_bioinformatic_pipeline.git
```

Dive to the newly created directory:
```
cd V4_SSU_ASV_bioinformatic_pipeline/
```

Adjust the number of CPU cores available on your machine (replace *xx* by the number of cores availalbe):
```
sed 's/^NCPUS=32$/NCPUS=xx/' script_V4.sh > my_script_V4.sh
```

Execute the example script:
```
bash ./my_script_V4.sh
```

Alternatively, open the files *my_script_V4.sh* and *script_dada2.sh* in your prefered text editor or in Rstudio, and execute steps one by one using copy & paste in a bash or R terminal.

Steps 1 to 6 and 15 to 17 are executed in a bash terminal.

Steps 7 to 14 are executed in R.
