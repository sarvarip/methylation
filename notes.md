- Symbolic link
    - ln -s /home/cmb-panasas2/sarvari ~/panfs
    - Ln -s /auto/cmb-05/qbio/sarvari

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

- Remove
    - rm sratoolkit.2.9.2-ubuntu64.tar.gz
    - Everything in subfolders: rm -rf sarvari$
    - rm -rf Python-3.4.3/

- Untar
    - Tar xzf sratoolkit.2.9.2-ubuntu64.tar.gz

- Request interactive computing job
    - Always use salloc for any job, no job should be done on the head node
    - salloc --partition=cmbr --mem=120GB --time=50:00:00 --ntasks=1 --cpus-per-task=16

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
    - export PATH=$PATH:~/methpipe/bin (add ~/methpipe/bin to the path)
    - export PATH=$PATH:~/sarvari/sratoolkit.2.9.2-ubuntu64/bin
    - export PATH=$PATH:~/miniconda3/bin
    - export PATH=$PATH:/usr/usc/matlab/R2018a/bin

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


- Sort the reads 
    - LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o mouse.mr.sorted mouse.mr
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
