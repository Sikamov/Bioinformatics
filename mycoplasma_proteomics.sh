#/storage/data2/micoplasma_proteomics/OV_706-735
#scp -P 2222 -i ssh /Users/kirillsikamov/Desktop/quantms_design.tsv sikamov@ssh.enteromics.com:/home/sikamov/micoplasma 
scp -P 2222 -i ssh /Users/kirillsikamov/Desktop/ATCC sikamov@ssh.enteromics.com:/storage/data2/sikamov/mycoplasma/



#NCBI 
#https://www.ncbi.nlm.nih.gov/genome/3075?genome_assembly_id=210048


#uniprot (не брал)
#https://www.uniprot.org/proteomes/UP000029711



###annotation

mkdir prefetch

prefetch SRR15616380
vdb-validate ./SRR15616380/ 
fasterq-dump SRR15616380/

mkdir fastqc
fastqc *.fastq -o fastqc

mkdir trim

#trimmomatic SE SRR15616380.fastq SRR15616380_2.fastq trim/pair1.fastq trim/unpair1.fastq trim/pair2.fastq trim/unpair2.fastq HEADCROP:3\ 
 TRAILING:20 LEADING:20 ILLUMINACLIP:trim/illumina_adapters.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:25

#fastqc trim/*.fastq -o fastqc

mkdir spades 

spades.py -s SRR15616380.fastq -o spades --careful --phred-offset 33



#-r GCF_009176875.1_ASM917687v1_genomic.fna.gz 
#/home/user12/genomics/tools/quast-5.2.0/

conda activate quastenv
mkdir quast_results
quast.py spades/scaffolds.fasta -o quast_results

scp -P 2222 -i ssh -r sikamov@ssh.enteromics.com:/storage/data2/micoplasma_proteomics/OV_706-735/references/quast_results/report.html /Users/kirillsikamov/Desktop/

conda activate bactaenv


bakta \
--db /storage/data2/sikamov/db \
--translation-table 4 \
-v --genus Mycoplasma \
--species "hominis" \
--strain H34 \
--locus-tag MHO \
--prefix Mh-H34_filt_bold \
--threads 8 \
--complete \
--compliant \
--proteins /storage/data2/sikamov/mycoplasma/GCF_000085865.1_ASM8586v1_genomic.faa \
--force \
--output /storage/data2/sikamov/mycoplasma/bakta_Mh-H34_filt_bold \
/storage/data2/sikamov/mycoplasma/genomes/Mh_H-34_filt_bold/assembly.fasta

bakta \
--db /disk/data/baktadb/db \
--translation-table 4 \
-v --genus Mycoplasma \
--species "hominis" \
--strain H34 \
--locus-tag MHO \
--prefix MHO \
--threads 10 \
--compliant \
--output $PWD/bakta_H34_contig \
--proteins $PWD/myco_atcc2314.faa \
$PWD/assembly.fasta



seqkit seq -a UP000002631_347256.fasta > UP000002631_347256.faa

bakta_proteins \
--db /storage/data2/sikamov/db \
--prefix Mh-H34_NCBI \
--output /storage/data2/sikamov/mycoplasma/bakta_Mh_H34_NCBI \
--proteins /storage/data2/sikamov/mycoplasma/UP000002631_347256.faa \
--threads 8 \
--force \
/storage/data2/sikamov/mycoplasma/GCF_000759385.1_ASM75938v1_protein.faa


############################

conda activate ngsenv

export JAVA_CMD="/usr/bin/java"
export JAVA_HOME="/usr/bin/java"
export TOWER_ACCESS_TOKEN="eyJ0aWQiOiA3NzA5fS43ZWNmMGE5ZjZkZjBhYTI2MTk2NzUyMjJjMGZlZDk0MTM0MWUzODE5"
export NXF_VER=23.04.0

#Mh_11 Mh_40

for i in  Mh_45 Mh_33 Mh_43 Mh_1862 Mh_12 Mh_X-37 Mh_7
do

mkdir ${i} 
cd ${i} 

nextflow run nf-core/quantms \
-profile docker \
--search_engines comet \
--search_engines msgf \
--local_input_type raw \
--root_folder /storage/data2/micoplasma_proteomics/OV_706-735 \
--input /storage/data2/micoplasma_proteomics/OV_706-735/mycoplasma_design_${i}.tsv \
--database /storage/data2/sikamov/mycoplasma/bakta_Mh-H34_filt_bold/Mh-H34_filt_bold.faa \
--add_decoys \
--max_precursor_charge 7 \
--enzyme Trypsin \
--fixed_mods 'Carbamidomethyl (C)' \
--variable_mods 'Oxidation (M)' \
--outdir /storage/data2/micoplasma_proteomics/OV_706-735/${i}  \
--max_cpus 8 \
--max_memory 30.GB \
-r 1.1.1 \
--labelling_type 'label free sample' \
--acquisition_method dda \
-with-tower \
--contrasts 2-1

telegram-send ${i}
cd ../
done


#scp -P 2222 -i ssh -r sikamov@ssh.enteromics.com:/storage/data2/micoplasma_proteomics/OV_706-735/Mh_*/msstats/mycoplasma_design_Mh_*  /Users/kirillsikamov/Desktop/МФТИ/Биоинформатика\ ФХМ/micoplasma/H34/ 

#scp -P 2222 -i ssh -r sikamov@ssh.enteromics.com:/storage/data2/micoplasma_proteomics/OV_706-735/Mh_*/proteomicslfq/*.csv /Users/kirillsikamov/Desktop/МФТИ/Биоинформатика\ ФХМ/micoplasma/H34/PCA/

#scp -P 2222 -i ssh -r sikamov@ssh.enteromics.com:/storage/data2/sikamov/mycoplasma/bakta_Mh_H34_filt_bold/Mh-H34_filt_bold.gff3 /Users/kirillsikamov/Desktop/


makeblastdb -in  UP000002631_347256.faa -out ATCC -dbtype prot -parse_seqids

blastn -query bakta_Mh-H34_PCM_ATCC23114/Mh-H34_PCM_ATCC23114.faa -db ATCC -outfmt "7 qacc sacc evalue qcovs pident" -out H34_vs_ATCC -perc_identity 95 -qcov_hsp_perc 90

../pgap.py -r -o pgap_MH-34_circle -g /storage/data2/sikamov/mycoplasma/genomes/Mh_H-34_cycle/Mh_H-34_cycle.fasta -s 'Mycoplasma hominis'