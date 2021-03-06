Project

Motility

- sed -i 's~,~\t~g' RMAdataMotility.csv
- grep "\S" RMAdataMotility.csv > Motility.tsv
- awk 'BEGIN{FS=OFS="\t"} {sub(/-.+/, "", $1)} 1' Motility.tsv > motility
- sed -i 's~"~~g' motility
- awk 'BEGIN{FS=OFS="\t"} {sub(/ .+/, "", $1)} 1' motility > motility.tsv
- vim edit to remove one weird line
- awk -F"\t" '{if (!($2=="" || $3=="" || $4=="")) print $0}' motility.tsv > motility.cleared
- awk -F"\t" '{if (NR>1) print $1"\t"100*$4/($2*$3); else print "Laboratory.ID\tMotility"}' motility.cleared > motility.tsv


General EWAS analysis

- library(data.table)
- location: /staging/as/andrewds/ewas/GSE55763_normalized_betas_clean.txt
- df <- fread(text = 'GSE55763_normalized_betas_clean.txt')

Aston smoking

- conda create -n methylation_r r bioconductor-minfi
- source activate methylation_r
- conda install r bioconductor-illuminahumanmethylation450kmanifest 
- dirs <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)[-1]
- rgSet <- read.metharray.exp(dirs)
- MSet.swan <- preprocessSWAN(rgSet)
- Mset.swan.betas <- getBeta(MSet.swan)

- choose 1 core in HPC if preprocessQuantile is run
- grSet <- preprocessQuantile(rgSet)
- grSet.beta <- getBeta(grSet)

- write.table(df, 'AgingNormalizedTogether.txt', sep='\t', quote=FALSE)

- multimethstat -v ./intervals/Human_Sperm_named.bed EPIC_hg19_probe_coords.bed ./input/AgingNormalizedTogether.txt -o AgingNorm_Human_Sperm_hmr.txt

Aging analysis

- data originally transferred to: sarvari@smithlab:/labdata/jenkins
- sarvari@smithlab.usc.edu

EPIC analysis

- Get promoters 1000 upstream and downstream of a gene 
- Use the Table browser for that: https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=728013001_Upd5TraWK5jaeJEmxkb0K01CKWaB
- concat them
- awk '$3 = $3+1000' proms.bed 
- awk '($6=="+") {print $1"\t"$2"\t"$3}' proms.bed | uniq > promstabs.bed
- sortBed -i promstabs.bed > prom.bed
- bedtools merge -i prom.bed > merged.bed
- sort -V -k 1,3 "merged.bed" -o proms_merged_sorted
- CpG islands: https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=578954849_wF1QP81SIHdfr8b0kmZUOcsZcHYr&clade=mammal&org=Human&db=hg38&hgta_group=regulation&hgta_track=knownGene&hgta_table=0&hgta_regionType=genome&position=chr9%3A133252000-133280861&hgta_outputType=primaryTable&hgta_outFileName=
- export PATH=$PATH:~/panfs/smithlab/methpipe/src/analysis
- awk 'BEGIN{k=0} {$4 = $4""k; k+=1; print}' cpgIslandExt_hg19_080913_good.bed > cpgIslandExt_hg19_080913_named.bed
- multimethstat -progress -v -o cpg_island_features_only.txt cpgIslandExt_hg19_080913_named.bed EPIC_hg19_probe_coords.bed methylation.txt
- Sorting according to multimethstats definition: LC_ALL=C sort -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o proms_merged_sorted_otherdef.bed proms_merged_sorted or simply sortBed -i promstabs.bed > prom.bed
- sort -V is to be used if bedtools merging is used after!! (otherwise bedtools error)
- awk 'BEGIN{k=0} {$4 = $1"_"$2"_"$3; print $0"\t0\t+"}' merged.bed > merged_sixcol.bed
- ~/panfs/methpipe_edit/methpipe/src/analysis/multimethstat -progress -v -o proms_features_only_fixed.txt merged_sixcol.bed EPIC_hg19_probe_coords.bed methylation_fixed.txt
- awk '{print $4}' merged_sixcol_sorted > colnames
- awk '{$1=$1+1; print $1}' gene_array > colidx %add one because in Python indexing starts from 0, but in R from 1
- df <- read.table("proms_features_only_fixed.txt", header = TRUE, sep = "\t", row.names = 1)
- To avoid shifting of columns, use instead: X <- read.table('file.txt', header=T)

Debug multimethstat

- grep -n chr6_72959032_72961032 merged_sixcol_sorted
- grep -n 6759 colidx
- head -n 6044 coord_lists_per_gene | tail -1
- awk '($4 == "chr6_72959032_72961032" || $4 == "chr9_115013537_115015537" || $4 == "chr15_30937316_30939316" || $4 == "chr15_90783140_90785140" || $4 == "chrY_9527708_9529708") {print $0}' merged_sixcol.bed > test.bed
- for i in `head -1 ../../multimethstat/methylation.txt `; do echo $i; done | grep cg114996
- head -3 methylation_fixed.txt | cut -d " " -f 1-5

Machine Learning prep

- awk '{ for(i=2;i<=10000;i++) {printf "%s ", $i} printf "\n"}' methylation.txt | awk '{for (i=1;i<=NF;i++) {if (NR==1) name[i] = $i; else total[i]+=$i; sq[i]+=$i*$i}} END {for(i=1;i<=NF;i++) {print name[i], sq[i]/NR-(total[i]/NR)**2}}'
- awk '{for (i=1;i<=NF;i++) {if (NR==1) name[i] = $i; else total[i]+=$i; sq[i]+=$i*$i}} END {for(i=1;i<=NF;i++) {print name[i], sq[i]/NR-(total[i]/NR)**2}}' methylation.txt > variances
- /usr/usc/R/3.5.0/bin/R
- readRDS('methylation.rds')
- variances <- apply(df,2,var)
- variances <- variances[,order(-variances)]
- df_red <- df[,names(variances[1:1000])]
- rownames(df_red) <- gsub("X", "", rownames(df_red))
- cut -f1 RMAdataAge.txt | cut -d "-" -f1 | cut -d "*" -f1 | cut -d "," -f1 > ids
- sed -i 's/"//g' ids
- cut -f2 RMAdataAge.txt > ages
- paste ids ages > agedat
- ages <- read.csv(file="agedat", header=TRUE, sep="\t")
- rownames(df_red) <- gsub("X", "", rownames(df_red))
- ages$Laboratory.ID <- NULL
- data <- merge(df_red, ages, by=0)
- rownames(data) <- data$Row.names
- data$Row.names <- NULL
- sed -i 's~n/a~~g' numsucc
- paste ids numeggs numsucc > temp
- awk -F"     " '{if (!($2=="" || $2==0 || NR==1)) print $3/$2; else print""}' temp > ratios
- sed -i 's~.*Endometrial.*~0~g' issue
- paste ratiodat issue > ratiodat_filtered
- awk -F"\t" '!($4 == 0) { print $1"\t"$2"\t"$3 }' ratiodat_filtered > filtered_ratiodat
- awk -F"\t" ' ($1=="cg02274263") {print NR}' sorted_EPIC_probe_coords
- awk -F"\t" ' ($1=="cg21078414") {print NR}' sorted_EPIC_probe_coords
- sed -n '629496,629746p' sorted_EPIC_probe_coords > SNORD_probe_coords
- awk {'print $1'} SNORD_probe_coords > SNORD115coords
- sed -i 's~.*#DIV/0.*~~g' malebmi.txt
- sed -i 's~0.*~~g' malebmidat
- awk '$2 > 40 && $5 == 0{print $0}' filtered_successdatratio | wc -l
- sed -i 's~,~\t~g' OldDataFemale_Age_Euploid_Rate.csv
- ~/panfs/ivfdat/IVF_regression/reproduce/final/euploid_analysis$ /usr/usc/R/3.5.0/bin/Rscript euploid_analysis.R features.txt metadata.txt methylation.txt Euploid.Rate EuploidAnalysis

Python codes / preparation / correlation validation

- cat 13_sites.bed | while read line; do echo "$line"; done | awk '{print $1"\011"$2}' >> 13_cols
- python extract.py -g 13_cols -p EPIC_hg19_probe_coords -m methylation.txt -o 13_reduced_methylation.txt (produces found_columns.txt!!!)
- python correlation.py -m1 aplicon_valid_df.txt -m2 13_reduced_methylation.txt -col found_cols.txt -o 13_correlation_full.txt
- python correlation.py -m1 well_b1_validation.txt -m2 13_reduced_methylation.txt -col found_cols.txt -o 13_correlation_b1.txt
- python correlation.py -m1 well_b2_validation.txt -m2 13_reduced_methylation.txt -col found_cols.txt -o 13_correlation_b2.txt
- python correlation.py -m1 well_b3_validation.txt -m2 13_reduced_methylation.txt -col found_cols.txt -o 13_correlation_b3.txt
- python sample_correlation.py -m1 aplicon_valid_df.txt -m2 13_reduced_methylation.txt -col found_cols.txt -o sample_correlation_full.txt
- python sample_correlation.py -m1 well_b3_validation.txt -m2 13_reduced_methylation.txt -col found_cols.txt -o sample_correlation_b3.txt
- python sample_correlation.py -m1 well_b2_validation.txt -m2 13_reduced_methylation.txt -col found_cols.txt -o sample_correlation_b2.txt
- python sample_correlation.py -m1 well_b1_validation.txt -m2 13_reduced_methylation.txt -col found_cols.txt -o sample_correlation_b1.txt
- cut -f2 correlation_full_doublecheck.txt 
- python locsort.py -i sites -i2 correlations -o sorted_sites -o2 sorted_correlations
- paste sorted_sites sorted_correlations > sorted_correlation_full_doublecheck.txt

- when getting the input from merge-methcounts, add in an extra tab at the beginning for Python pandas to read the file well, also need to edit the weird row and column names (also follow other steps in Untitled.ipynb)
- python correlation.py -m1 aplicon_valid_df.txt -m2 13_reduced_methylation.txt -col found_cols.txt -o 13_correlation_full.txt -d1 amplicon_13_SAS.csv -d2 epic_13_SAS.csv

Amplicon validation

- ls *.fastq > filenames.txt
- ls *[0-9].fastq > batch1names.txt
- for i in *.fq; do mv "$i" "$(basename "$i" .fq).fastq"; done
- mv batch1names.txt /staging/as/sarvari/ivfdat
- ~/panfs/bismark/Bismark-0.22.1/bismark_genome_preparation --parallel 32  --verbose chroms
- xargs --arg-file=command_list0.txt --max-procs=32 --replace --verbose /bin/sh -c "{}"
- xargs --arg-file=command_list_p1.txt --max-procs=13 --replace --verbose /bin/sh -c "{}"
- xargs --arg-file=command_list_p2.txt --max-procs=13 --replace --verbose /bin/sh -c "{}"
- xargs --arg-file=command_list_onehalf.txt --max-procs=15 --replace --verbose /bin/sh -c "{}"
- xargs --arg-file=command_list2.txt --max-procs=32 --replace --verbose /bin/sh -c "{}"
- xargs --arg-file=command_list3.txt --max-procs=32 --replace --verbose /bin/sh -c "{}"
- xargs --arg-file=command_list4.txt --max-procs=32 --replace --verbose /bin/sh -c "{}"
- Much better way: ls *_cutadapt.mr.sorted | xargs -P 16 -n 1 -I % /home/rcf-47/andrewds/as/code/for_peter/methpipe/bin/duplicate-remover   (uses the 16 cores for one job though, automatically named output)
- To use multiple salloc cluster jobs:
- for i in *.mr.sorted; do echo $(basename $i ".mr.sorted").meth; done | xargs snakemake -p -s methcounts.smk -j 200 --cluster "{params.grid_opts}" --rerun-incomplete --latency-wait 60
- for i in *.fastq; do echo $(basename $i ".fastq")_trimmed.fq; done | xargs snakemake -p -s cutadapt.smk -j 20 --cluster "{params.grid_opts}" --rerun-incomplete --latency-wait 60
- for i in *.meth; do echo $(basename $i ".meth").meth.highcoverage; done | xargs snakemake --use-singularity -p -s methfilter.smk -j 20 --cluster "{params.grid_opts}" --rerun-incomplete --latency-wait 60
- awk -F "\t" '$6 > 1000' 61578_S65_cutadapt.meth > 61578_S65_cutadapt.meth.highcoverage
- awk '$1 == "chr10" && $2 == "45641678" {print $5,$6}' 59160_S209_cutadapt.meth.highcoverage
- for i in *.highcoverage; do wc -l "$i" >> linecounts.txt; done
- awk '$1 == 0 {print $2}' linecounts.txt 
- for i in *.highcoverage; do awk '$1 == "chr10" && $2 == "45641678" {print $5,$6, FILENAME}'  "$i" >> summary_s1.txt; done 
- /home/rcf-47/andrewds/as/code/for_peter/methpipe/bin/selectsites -v 13_sites.bed 57525.meth
- for i in *.highcoverage; do awk '$6 > 65536 {print FILENAME}' "$i" >> biggerthan65.txt; done
- for i in *.meth; do echo $(basename $i ".meth").meth.roi; done | xargs snakemake -p -s methroi.smk -j 20 --cluster "{params.grid_opts}" --rerun-incomplete --latency-wait 60
- for i in *.mr.sorted; do echo $(basename $i ".mr.sorted").bsrate; done | xargs snakemake -p -s bsrate.smk -j 20 --cluster "{params.grid_opts}" --rerun-incomplete --latency-wait 60
- for i in *.meth; do echo $(basename $i ".meth").levels; done | xargs snakemake -p -s levels.smk -j 10 --cluster "{params.grid_opts}" --rerun-incomplete --latency-wait 60
- awk '$1 > 57235792 {print$1,$2}' methlinecounts.txt > aberrantmethfiles.txt
- cat filenametruncmeth.txt | while read line; do echo "$line"; done | xargs snakemake -p -s methcounts.smk -j 20 --cluster "{params.grid_opts}" --rerun-incomplete --latency-wait 60 (don't do this snakemake automatically knows which ones to rerun)
- cat filenametruncmeth.txt | while read line; do rm -rf "$line"; done
-  /home/rcf-47/andrewds/as/code/for_peter/methpipe/bin/merge-methcounts -v -h -o amplicon_validation_mergedmeths_nodasht.txt *.roi
- for i in *.sam; do t=$(basename $i “.sam”); n_mapped=$(grep -v ^@ ${t}.sam | cut -f 1 | uniq | wc -l); n_seqd=$(cat ${t}.fastq | wc -l); echo $t $n_mapped $n_seqd; done | awk ‘{print $1,$2,$3/4}’ > seqd_and_mapped.txt
- /home/rcf-47/andrewds/as/code/for_peter/methpipe/bin/merge-methcounts -t *.roi > amplicon_validation_mergedmeths.txt
- cat batch2names.txt | while read line; do echo $(basename $line ".fastq")_cutadapt.fastq; done >> batch2namescutadapt.txt 
- for i in *.fastq; do awk 'NR%4==2 {print length}' "$i" >> merged_seqlen.txt; done
- for i in `cat //staging/as/andrewds/rma/batch_3.txt | cut -f 1`; do echo $i `grep -v ^@ ${i}*sam | wc -l` `cat ${i}*mr | wc -l`; done | awk '{print $1,$2,$3,$2==$3}' >> ~sarvari/panfs/ivfdat/good_tomr_conversions.txt
- for i in `grep 0$  ~sarvari/panfs/ivfdat/good_tomr_conversions.txt | cut -f 1`; do find ~sarvari/staging/ -maxdepth 1 -name ${i}\*_bismark_bt2.bam; done > redo_tomr.txt
- grep 0$ good_tomr_conversions.txt | awk '{print $1}'
- for i in `grep 0$  ~sarvari/panfs/ivfdat/good_tomr_conversions.txt | cut -f 1`; do find ~sarvari/staging/ -maxdepth 1 -name ${i}\*.mr; done > redo_sort.txt
- for i in *.bsrate;do awk 'NR==1 {print FILENAME,$5}' "$i" >> bsrate_summary.txt; done
- sed -e "s/_S[0-9]\+_cutadapt.bsrate//g" ~sarvari/staging/bsrate/bssum.txt
- sed -e "s/.bsrate//g" ~sarvari/staging/bsrate/bssum.txt >> summary.txt
- sort -k1 -n summary.txt >> sorted_summary.txt

- awk 'NR%4==2 { thislen=length($0); totlen+=thislen}
END { print FILENAME,  4*totlen/NR } ' 60098_S186_cutadapt.fastq

- for i in *_S*_cutadapt.fastq; do awk 'NR%4==2 { thislen=length($0); totlen+=thislen}
END { print FILENAME,  4*totlen/NR } ' "$i"; done >> fraglen

- for i in *[0-9].fastq; do awk 'NR%4==2 { thislen=length($0); totlen+=thislen}
END { print FILENAME,  4*totlen/NR } ' "$i"; done >> fraglen

- cat fraglen1 fraglen23 > fraglen

General

- find /home/cmb-panasas2/sarvari -type d -exec chmod 775 \{} \; to change permissions for all files in panfs such that they can be edited by group members
- R install.packages
    - move R directory at login node (where packages are automatically installed) to panasas and create a symbolic link from panasas to the login directory

- Conda
    - conda installation is taking long if only one environment is used
    - conda create -n methpipe 
    - conda activate methpipe
    - conda install -c psi4 gcc-5 
    - conda install -c conda-forge gsl 
    - conda deactivate

- Jupyter notebook
    - jupyter notebook --no-browser --port=8889
    - ssh -N -f -L localhost:8888:localhost:8889 sarvari@hpc-cmb.usc.edu
    - open web browser with URL: localhost:8888
    - copy token 
    - If cannot listen error, kill all processes on 8888: kill -9 $(lsof -t -i:8888)
    
- Hardware info
    - cat /proc/meminfo
    
- C++ compiler location
    - export PATH=$PATH:/usr/usc/gnu/gcc/5.3.0/bin
    - Still does not work because the path to the earlier version 
    - export PATH="/usr/usc/gnu/gcc/5.3.0/bin:/home/rcf-40/sarvari/panfs/miniconda3/bin:/home/rcf-40/sarvari/panfs/miniconda3/condabin:/usr/lib64/qt-3.3/bin:/us/bin:/usr/bin:/usr/local/sbin:/usr/sbin"
    - set it to be path in .bashrc (just copy it to the end of file) vim ~/.bashrc
    
- Misc
    - Jump to front of end of line: CTRL-A and CTRL-E
    - Reverse search of commands: CTRL-C and then CTRL-R
    - Paste: Middle mouse button
    - Put running process in background: CTRL-Z and then type bg
    - Type fg to retrieve it
    - Check running jobs: squeue -u sarvari
    
- Bash commands for all files in a directory
    - for i in *.bz2; do ln -s ~/panfs/ivfdat/amplicon_validation/20190415/fastq/"$i" /staging/as/sarvari/ivfdat; done
    - for i in *.bz2; do bzip2 -d "$i"; done 
    - for i in *.fq; do mv "$i" "$(basename "$i" .fq).fastq"; done
    - Going through lines of a file:  IFS=$'\n'; for j in $(cat ./command_list.txt); do echo "$j"; done

- Symbolic link
    - ln -s /home/cmb-panasas2/sarvari ~/panfs
    - cat batch3names.txt | while read line; do ln -s ~/panfs/ivfdat/amplicon_validation/20190424/fastq/"$line" ~/staging/fastq; done;
    - IFS=$'\n'; set -f; for i in $(cat < "mrsortedhead.txt"); do ln -s ~/staging/"$i" ~/staging/snakemake_test; done;
    
- Text editor
    - Use https://github.com/spf13/spf13-vim
    
- Copy a specific column into a new file
    - awk {'print $2'} EPIC_hg19_probe_coords > EPIC_chroms

- Count lines
    - wc -l mouse.mr.dremove

- Permission for a library (data)
    - Chmod 777 data

- SCP (-r for folder)
    - Scp -r mengzhou sarvari@hpc-cmb.usc.edu:~/panfs/

- Download
    - Copy link address from website
    - Use wget and paste after
    - For SRA download use this database: ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR306/SRR306421/
    - Human genome: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/
    - Use rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ . to download

- Remove
    - rm sratoolkit.2.9.2-ubuntu64.tar.gz
    - Everything in subfolders: rm -rf sarvari$
    - rm -rf Python-3.4.3/

- Untar/Unzip
    - Tar xzf sratoolkit.2.9.2-ubuntu64.tar.gz
    - gzip -d well_covered_merged.tsv.gz

- Request interactive computing job
    - Always use salloc for any job, no job should be done on the head node
    - salloc --partition=cmb --mem=120GB --time=50:00:00 --ntasks=1 --cpus-per-task=16 --account=lc_as

- Open big file
    - head mouse.mr > mouse.mr.head
    - vim mouse.mr.head
    
- GIT
    - git clone https://github.com/sarvarip/methylation.git
    - git add notes.md
    - git commit -am "Added notes.md"
    - git push
    - git pull

- VIM
    - I: insert
    - esc: command mode
    - :w save
    - :q! exit
    - :wq to exit and save
    - :set paste from CTRL SHIFT V pasting (so that it does not screw up with the autoindent)
    - :set nopaste to change it back
    - :vsp filename to vertically split window and can easily copy and paste using 'y' and 'p'

- Emacs
    - Ctrl-x and CTRL-s to save
    - Ctrl-x and CTRL-c to exit
    - Press Q if something weird comes up
    - New script: emacs script1.sh
    - Run: sbatch script1.sh

- Move
    - mv SRR5381655_* ~/sarvari/ (star means move all with such beginning)

- Changing the environmental path (MUST DO IT EVERY TIME)
    - echo $PATH (let’s you see the current environmental variables)
    - export PATH=$PATH:~/panfs/smithlab/methpipe/src/analysis/
    - export PATH=$PATH:~/sarvari/sratoolkit.2.9.2-ubuntu64/bin
    - export PATH=$PATH:~/miniconda3/bin
    - export PATH=$PATH:/usr/usc/matlab/R2018a/bin
    - export PATH=$PATH:/usr/usc/R/3.5.0/bin

- Create new environmental variable
    - export SMITHLAB_CPP=~/panfs/smithlab/smithlab_cpp

- Check sizes and read files
    - Ls -l (check size)
    - Less (read file, exit with Q)

- Check status of job and cancel it
    - squeue -u sarvari (username)
    - scancel jobid

- Put current job in background
    - CTRL-Z and then type bg
    - To put it back in foreground type fg

- TMUX
    - Tmux new -s name to create new session
    - CTRL-B C to create new window within session
    - Then to change windows, use CTRL-B and the number of the window you wanna jump to 
    - CTRL-B d to detach the session
    - Tmux a -t name to attach the session
    - Exit to close down the session
    - Scroll: CTRL-B and Page Up button

- Install bioconda
    - After installing conda: 
    - https://bioconda.github.io/
    - conda install bowtie2


- Using matlab
    - /usr/usc/matlab/R2018b/bin/matlab


- Login from Ubuntu
    - Ssh sarvari@hpc-cmb.usc.edu
    
Pipeline

- Methpipe installation
    - git clone --recursive https://github.com/smithlabcode/methpipe.git
    - git clone --recursive https://github.com/smithlabcode/smithlab_cpp.git
    - use g++ version 5+ (see above how to get it)
    - rename htslib_wrapper.*pp such that it does not have any extension
    - make -f original_makefile.mk static in smithlab_cpp directory
    - create new environmental variable, SMITHLAB_CPP (see above)
    - make -f original_makefile.mk in mehpipe directory 
    
    - (1) Mac users will have to accept the Xcode license before using any of the command line tools related to development.
    - (2) GSL must be installed
    - (3) Brew must be installed (needed for GSL)
    - (4) The default location (and the right one) for brew is in /usr/local, but the user does not have write access    (solution: recursively change ownership of /usr/local/* to the user)
    - (5) The default location to install brew is /usr/local/homebrew, which required that we create symbolic links (using ln). This should be changed to just /usr/local instead.
    - (6) We need to make sure that by default the htslib_wrapper files are excluded so users don’t need to worry about those.
    - Also, it would be really great to find out why g++ 4.2.1 worked for Tim, but not others in the Smith lab.
    
- Snakemake file for the whole pipeline
    - /home/rcf-47/andrewds/as/code/MethBase/MethBaseConfig

- Downloads
    - Download the genome in a format: chromFa.tar.gz
    - Untar it: tar -xvzf chromFa.tar.gz


- Makedb
    - ~/walt/bin/makedb -c mousechroms -o mouse.dbindex


- Fastq-dump 
    - Location: sarvari/sratoolkit.2.9.2-ubuntu64/bin
    - Command: ./fastq-dump –split-files ../../SRR5381655.sra


- Cutadapt
    - /home/rcf-40/sarvari/.local/bin/cutadapt --help
    - Adaptor sequence: AGATCGGAAGAGC
    - cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o 1M-2wks_1_cutadapt.fastq -p 1M-2wks_2_cutadapt.fastq 1M-2wks_1.fastq 1M-2wks_2.fastq


- Paired-end read mapping with walt
    - sarvari@hpc3339:/panfs/cmb-panasas2/sarvari/data/Mouse_Kidney/reads$ 
    - ~/walt/bin/walt -t 16 -i mouse.dbindex -1 1M-2wks_1_cutadapt.fastq -2 1M-2wks_2_cutadapt.fastq -o mouse.mr

- Download Mouse full genome in one file
    - UCSC genome browser mm10.2bit
    - Download TwoBitToFa from UCSC Genome Browser
    - Make it executable by Chmod 777  TwoBitToFa
    - Then ~/panfs/twoBitToFa mm10.2bit mousegenome.fa


- Clone and checkpoint to get walt2
    - git clone --recursive git@github.com:andrewdavidsmith/walt2.git
    - cd walt2
    - git checkout bf8046a5bd6f30cd7c8a8b4f507a94bcd48520fb
    - request computing job!!
    - ~/panfs/walt2/bin/walt2idx -v mousegenome.fa mouse.walttwoindex
    - ~/panfs/walt2/bin/walt2 -t 16 -i mouse.walttwoindex -o mouse_walt2.mr 1M-2wks_1_cutadapt.fastq 1M-2wks_2_cutadapt.fastq



- Bismark
    - ~/panfs/bismark/Bismark-0.22.1/bismark_genome_preparation --verbose mousechroms
    - ~/panfs/bismark/Bismark-0.22.1/bismark --multicore 16 --genome /panfs/cmb-panasas2/sarvari/mousechroms -1 1M-2wks_1_cutadapt.fastq -2 1M-2wks_2_cutadapt.fastq
    - to-mr -m bismark -v -o mouse_bismark.mr 1M-2wks_1_cutadapt_bismark_bt2_pe.bam


- Sort the reads 
    - LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -S120G --parallel=16 -o mouse_parallel.mr.sorted mouse.mr
    - LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -S120G --parallel=16 -o mouse_walt2.mr.sorted mouse_walt2.mr

- Remove duplicates
    - Duplicate remover is in the methpipe folder
    - export PATH=$PATH:~/methpipe/bin
    - duplicate-remover -S mouse_dremove_stat.txt -o mouse.mr.dremove mouse.mr.sorted
    - duplicate-remover -S mouse_parallel_dremove_stat.txt -o mouse_parallel.mr.dremove mouse_parallel.mr.sorted
    - duplicate-remover -S mouse_walt2_dremove_stat.txt -o mouse_walt2.mr.dremove mouse_walt2.mr.sorted

- Methcounts
    - methcounts -c ~/panfs/mousechroms -o mouse.meth mouse.mr.dremove
    - methcounts -c ~/panfs/mousechroms -o mouse_walt2.meth mouse_walt2.mr.dremove

- Check bsrate and levels and see if they make sense
    - bsrate -c mousechroms -o mouse.bsrate mouse.mr.dremove
    - bsrate -c /panfs/cmb-panasas2/sarvari/mousechroms -o mouse_walt2.bsrate mouse_walt2.mr.dremove
    - levels -o mouse.levels mouse.meth
    - levels -o mouse_walt2.levels mouse_walt2.meth


- Extracting and merging symmetric CpG methylation levels
    - symmetric-cpgs -o mouse_CPG.meth mouse.meth
    - symmetric-cpgs -o mouse_walt2_CPG.meth mouse_walt2.meth

- Hypomethylated regions
    - hmr -o mouse.hmr mouse_CPG.meth
    - hmr -o mouse_walt2.hmr mouse_walt2_CPG.meth

- Epiread format
    - methstates -c chroms -o SRR5381655_1.epiread SRR5381655_1_mapping.mr.dremove
    - methstates -c /panfs/cmb-panasas2/sarvari/mousechroms -o mouse_walt2.epiread mouse_walt2.mr.dremove

- Allelically methylated regions 
    - amrfinder -o SRR5381655_1.amr -c chroms SRR5381655_1.epiread
    - amrfinder -o mouse_walt2.amr -c /panfs/cmb-panasas2/sarvari/mousechroms mouse_walt2.epiread
    
- Roimethstat
    - LC_ALL=C sort -k 1,1 -k 3,3n -k 2,2n -k 6,6 -S120G --parallel=16 -o mouseproms.bed.sorted mouseproms.bed
    - roimethstat -o mouse_walt2.methstat ~/panfs/mouseproms.bed.sorted mouse_walt2_CPG.meth
    - use -L in roimethstat if there is a lot of memory requested
    
Visualization

    - Downloads:
        - wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig
        - wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v369/fetchChromSizes
        - wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v369/bedToBigBed
    - Permissions:
        - chmod 755 bedToBigBed
        - chmod 755 fetchChromSizes
        - chmod 755 wigToBigWig
    - Fetch chromosome sizes
        - ./fetchChromSizes hg38 > hg38.chrom.sizes
    - Convert the methcounts file to a bed format
        - awk -v OFS=”\t” ‘{print $1, $2, $2+1, $4”:”$6, $5, $3}’  SRR5381655_1_CPG.meth > SRR5381655_1_CPG.meth.bed
        - Rewrite all ‘’ and “” because Word formats them automatically in the wrong way
    - Create a bw track from the bed format methcounts output
        - cut -f 1-3,5 SRR5381655_1_CPG.meth.bed | ./wigToBigWig  /dev/stdin hg38.chrom.sizes SRR5381655_1_CPG.meth.bw
    - Make BigBed tracks for HMRs and AMRs
        - awk -v OFS=”\t” ‘{print $1,$2,$3,$4,int($5)}’ SRR5381655_1_CPG.hmr > SRR5381655_1_CPG.hmr.rounded
        - Rewrite all ‘’ and “” because Word formats them automatically in the wrong way
        - ./bedToBigBed SRR5381655_1_CPG.hmr.rounded hg38.chrom.sizes SRR5381655_1_CPG.hmr.bb
    - Visualize:
        - Make static IP and insert file links, e.g. http://128.125.86.68/tracks/setup.txt
        - Go to Genome Browser, Custom tracks, insert link and click go
        - Set Methylation to Full and HMR to dense
        - Right click on methylation data and configure methylation such that the windowing function is only mean (and not mean and whiskers)
        - To get the link, need to have a Desktop PC, or ask someone with one to read the files from your folder on the cluster, but then you need to make your files readable
        - Go to home folder and write chmod 755 sarvari
        - Then cd into your folder and do chmod 755 sarvari 
    - Result: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A174269623-248956422&hgsid=704058207_dJOlRmATUzWVAYFLojGxRsNghcww


Data:
    - The data I used for the above analysis was https://www.ncbi.nlm.nih.gov/sra?term=SRX2676780 which turned out to be RRBS (Reduced Representation Bisulfite Sequencing) data instead of WGBS, so the pipeline is not applicable

    - Get correct data: e.g. sperm data from AD Smith paper: https://www.sciencedirect.com/science/article/pii/S0092867411009421 Then see the Accession Number, which is GSE30340 and then look for this in the gene expression omnibus: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM752295 and then use both SRA samples and merge them after using Walt on all of them and continue with the pipeline

https://www.godaddy.com/garage/how-to-install-and-configure-python-on-a-hosted-server/
Sbatch script:

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mem=40000M
#SBATCH --cpus-per-task=16
#SBATCH -J test40G
#SBATCH -e home/rcf-40/sarvari/errorlog.txt
#SBATCH -o home/rcf-40/sarvari/screenoutput.txt
#SBATCH -p cmbr

cd /home/rcf-40/sarvari/sarvari
~/walt/bin/walt -t 16 -i hg38.dbindex -r SRR5381655_1_cutadapt.fastq -o SRR5381\
655_test40G__mapping.mr
