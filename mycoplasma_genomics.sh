
#/mnt/yandex.disk/Downloads/bacteria_myco_ecoli_gertsen_run_15-07-23/mhominis
conda activate fastqcenv
fastqc data/*.fq.gz -o qc/
fastqc data/Nanopore_H34_data/*.fastq.gz -o qc/

conda activate ngsenv
cd qc/
multiqc .
telegram-send --file multiqc_report.html 

cat data/Nanopore_H34_data/*.fastq.gz > data/Nanopore_H-34.fastq.gz

#Mh_7 Mh_43 Mh_11 Mh_40 Mh_33 Mh_12 Mh_1862 Mh_45


#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ngsenv


unicycler -1 data/V350171692_Run204_L02_36_1.fq.gz \
    -2 data/V350171692_Run204_L02_36_2.fq.gz \
    -l data/Nanopore_H-34_filt.fastq.gz \
    --mode bold \
    -o Mh_H-34_filt_bold/

telegram-send unicycler_H-34


for i in 28 29 30 114 32 33 34 35 
do

mkdir ${i}

unicycler -1 data/V350171692_Run204_L02_${i}_1.fq.gz \
    -2 data/V350171692_Run204_L02_${i}_2.fq.gz \
    -o ${i}/
telegram-send unicycler_${i}

done






for i in 28 29 30 114 32 33 34 35 
do

cd ${i}/bakta_${i}/
telegram-send --file ${i}.tsv
cd ../../


done

#scp -P 2222 -i ssh -r sikamov@ssh.enteromics.com:/storage/data2/sikamov/mycoplasma/genomes/Mh_H-34/*.gfa /Users/kirillsikamov/Desktop/


filtlong -1 V350171692_Run204_L02_36_2.fq.gz -2 V350171692_Run204_L02_36_2.fq.gz \
    --min_length 1000 --keep_percent 90 --target_bases 600000000 \
    --trim --split 100 \
    --mean_q_weight 10 Nanopore_H-34.fastq.gz | gzip > Nanopore_H-34_filt.fastq.gz


scp -P 2222 -i ssh -r sikamov@ssh.enteromics.com:/storage/data2/sikamov/mycoplasma/genomes/Mh_H-34_cycle/sorted_* /Users/kirillsikamov/Desktop/Mh_H-34_cycle/


# Индексация референтного генома для BWA и выравнивание коротких ридов
bwa index Mh_H-34_cycle.fasta
bwa mem Mh_H-34_cycle.fasta ../data/V350171692_Run204_L02_36_1.fq.gz ../data/V350171692_Run204_L02_36_2.fq.gz > short_alignment.sam

# Индексация референтного генома для Minimap2 и выравнивание длинных ридов
minimap2 -d Mh_H-34_cycle.mmi Mh_H-34_cycle.fasta
minimap2 -a Mh_H-34_cycle.mmi ../data/Nanopore_H-34_filt.fastq.gz > long_alignment.sam

# Преобразование коротких ридов SAM в формат BAM, их сортировка и индексация
samtools view -S -b short_alignment.sam > short_alignment.bam
samtools sort short_alignment.bam -o sorted_short_alignment.bam
samtools index sorted_short_alignment.bam

# Преобразование длинных ридов SAM в формат BAM, их сортировка и индексация
samtools view -S -b long_alignment.sam > long_alignment.bam
samtools sort long_alignment.bam -o sorted_long_alignment.bam
samtools index sorted_long_alignment.bam

echo "Alignment complete. Please load the following files into IGV:"
echo "Reference: Mh_H-34_cycle.fasta"
echo "Short reads alignment: sorted_short_alignment.bam"
echo "Long reads alignment: sorted_long_alignment.bam"




#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate canuenv

# Путь к данным
LONG_READS="data/Nanopore_H-34_filt.fastq.gz"
SHORT_READS_1="data/V350171692_Run204_L02_36_1.fq.gz"
SHORT_READS_2="data/V350171692_Run204_L02_36_2.fq.gz"

# Предварительная сборка с использованием длинных ридов
canu -p m_hominis -d canu_output genomeSize=0.7m -nanopore $LONG_READS

conda activate trycyclerenv
# Сборка с использованием Trycycler
trycycler cluster --reads $LONG_READS --assemblies canu_output/m_hominis.contigs.fasta --out_dir trycycler_output
trycycler reconcile --reads $LONG_READS --cluster_dir trycycler_output/cluster01 --out_dir trycycler_output/reconcile01

conda activate ngsenv
# Отображение и очистка сборки с использованием коротких ридов
minimap2 -ax sr trycycler_output/reconcile01/consensus.fasta $SHORT_READS_1 $SHORT_READS_2 > mapped.sam
racon $SHORT_READS_1 $SHORT_READS_2 mapped.sam trycycler_output/reconcile01/consensus.fasta > racon_output.fasta

# Использование Pilon для дополнительной коррекции сборки
pilon --genome racon_output.fasta --frags $SHORT_READS_1 --frags $SHORT_READS_2 --output pilon_output

# Оценка качества сборки
quast pilon_output.fasta -o quast_report

# Завершение скрипта
telegram-send готово!


###################################################

cd /Users/kirillsikamov/Desktop/genomes/
for i in 28 29 30 114 32 33 34 35 36 
do
    mkdir ${i}
    scp -i /Users/kirillsikamov/ssh sikamov.kv@calc.cod.phystech.edu:/home/common/sikamov.kv/mycoplasma/genomes/${i}/assembly.fasta /Users/kirillsikamov/Desktop/genomes/${i} 
    cd /Users/kirillsikamov/Desktop/genomes/
    #Sikamov2811
done




for i in 28 29 30 114 32 33 34 35 36 
do
    cd /storage/data2/sikamov/mycoplasma/genomes/MIPT/genomes/${i}
    /storage/data2/sikamov/pgap.py \
        -r -o pgap_${i} \
        -g assembly.fasta \
        -s 'Mycoplasma hominis' \
        --cpus 8
    telegram-send ${i}
    cd /storage/data2/sikamov/mycoplasma/genomes/MIPT/genomes/
done

cd /storage/data2/sikamov/mycoplasma/genomes/MIPT/genomes/36
/storage/data2/sikamov/pgap.py \
        -r -o pgap_36 \
        -g assembly.fasta \
        -s 'Mycoplasma hominis' \
        --cpus 8



###############################################################Unicycler

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=RT  # Замените YOUR_PARTITION на нужное значение
#SBATCH --job-name=unicycler_run
#SBATCH --mail-user=sikamov.kv@phystech.edu  # Замените на вашу почту
#SBATCH --comment="Unicycler assembly"

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate ngsenv

cd ~/mycoplasma/genomes/

mkdir 36
unicycler -1 ~/mycoplasma/data/V350171692_Run204_L02_36_1.fq.gz \
    -2 ~/mycoplasma/data/V350171692_Run204_L02_36_2.fq.gz \
    -l ~/mycoplasma/data/Nanopore_H-34_filt.fastq.gz \
    --mode bold \
    -o 36

for i in 28 29 30 114 32 33 34 35
do
    mkdir ${i}
    unicycler -1 ~/mycoplasma/data/V350171692_Run204_L02_${i}_1.fq.gz \
              -2 ~/mycoplasma/data/V350171692_Run204_L02_${i}_2.fq.gz \
              -o ${i}
done



############################################Prodigal


#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=RT  # Замените YOUR_PARTITION на нужное значение
#SBATCH --job-name=prodigal
#SBATCH --mail-user=sikamov.kv@phystech.edu  # Замените на вашу почту
#SBATCH --comment="Prodigal"

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate ngsenv
cd ~/mycoplasma/genomes/
for i in 28 29 30 114 32 33 34 35 36 
do
    cd ${i}
    mkdir prodigal_${i}
    prodigal -i assembly.fasta -a prodigal_${i}/proteins_${i}.faa -o prodigal_${i}/genes_${i}.gbk -p meta 
    cd ~/mycoplasma/genomes/
done

#!/bin/bash

# Путь к основной папке
base_dir=~/mycoplasma/genomes

# Список образцов
samples=(28 29 30 114 32 33 34 35 36)  # Добавьте все ваши образцы здесь

for sample in "${samples[@]}"; do
    file_path="$base_dir/$sample/prodigal_$sample/proteins_$sample.faa"
    temp_file="$base_dir/$sample/prodigal_$sample/temp.faa"

    # Добавляем название образца к идентификатору каждого белка
    sed "s/>/>${sample}_/" $file_path > $temp_file

    # Заменяем исходный файл модифицированным
    mv $temp_file $file_path
done


wc -l ~/mycoplasma/genomes/*/prodigal_*/proteins*

cat ~/mycoplasma/genomes/*/prodigal_*/proteins_*.faa > all_proteins.faa

srun --partition=RT --ntasks=1 --cpus-per-task=4 cd-hit -i all_proteins.faa -o clustered_proteins.faa -c 0.99 -n 5

# 1. Отфильтровываем кластеры
awk '/>Cluster/{if(f && c==9)print s; s=$0;f=0;c=0;next} /28|29|30|114|32|33|34|35|36/{f=1;c++} {s=s"\n"$0}' clustered_proteins.faa.clstr > filtered_clusters.faa.clstr

# 2. Извлекаем идентификаторы белков из отфильтрованных кластеров
grep -oP '>(\w+)' filtered_clusters.faa.clstr | sed 's/>//' > selected_proteins.txt

# 3. Оставляем только те белки, которые соответствуют выбранным идентификаторам
grep -A 1 -Fwf selected_proteins.txt clustered_proteins.faa > selected_proteins.faa

#######################################################PGAP

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=RT  # Замените YOUR_PARTITION на нужное значение
#SBATCH --job-name=PGAP
#SBATCH --mail-user=sikamov.kv@phystech.edu  # Замените на вашу почту
#SBATCH --comment="PGAP"

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate ngsenv
cd ~/mycoplasma/genomes/
for i in 28 29 30 114 32 33 34 35 36 
do
    cd ${i}
    ~/pgap/pgap.py \
        -r -o pgap_${i} \
        -g assembly.fasta \
        -s 'Mycoplasma hominis'
    cd ~/mycoplasma/genomes/
done


#rm -rf ~/mycoplasma/genomes/*/pgap_*


srun --partition=RT --ntasks=1 --cpus-per-task=8 ~/pgap/pgap.py \ -r -o pgap_36 \ -g assembly.fasta \ -s 'Mycoplasma hominis'


######################################################################BAKTA

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=RT  # Замените YOUR_PARTITION на нужное значение
#SBATCH --job-name=bakta
#SBATCH --mail-user=sikamov.kv@phystech.edu  # Замените на вашу почту
#SBATCH --comment="bakta"

source ~/miniconda3/etc/profile.d/conda.sh 

declare -A map
map[28]="7"
map[29]="43"
map[30]="11"
map[114]="40"
map[32]="33"
map[33]="12"
map[34]="1862"
map[35]="45"
map[36]="H-34"

for i in 28 29 30 114 32 33 34 35 36 
do
    j=${map[$i]}

    conda activate baktaenv
    bakta \
    --db ~/baktadb/db \
    --translation-table 4 \
    -v --genus Mycoplasma \
    --species "hominis" \
    --locus-tag MHO_${j}\
    --prefix MHO_${j} \
    --threads 8 \
    --complete \
    --compliant \
    --proteins ~mycoplasma/references/GCF_000085865.1_ASM8586v1_genomic.faa \
    --force \
    --output ~/mycoplasma/genomes/${i}/bakta_${i} \
    ~/mycoplasma/genomes/${i}/assembly.fasta

    conda activate ngsenv
    telegram-send bakta_${j}

done


conda activate baktaenv
    bakta \
    --db ~/bakta/db \
    --translation-table 4 \
    -v --genus Mycoplasma \
    --species "hominis" \
    --locus-tag MHOH34\
    --prefix MHO_H-34 \
    --threads 4 \
    --complete \
    --compliant \
    --proteins ~mycoplasma/references/GCF_000085865.1_ASM8586v1_genomic.faa \
    --force \
    --output ~/mycoplasma/genomes/genomes/36/bakta_36 \
    ~/mycoplasma/genomes/36/assembly.fasta




#######################################dfast

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=RT  # Замените YOUR_PARTITION на нужное значение
#SBATCH --job-name=dfast
#SBATCH --mail-user=sikamov.kv@phystech.edu  # Замените на вашу почту
#SBATCH --comment="dfast"
#SBATCH -o ~/mycoplasma/slurm_out/slurm/%j.out
#SBATCH -e ~/mycoplasma/slurm_out/slurm/%j.err

source ~/miniconda3/etc/profile.d/conda.sh 

declare -A map
map[28]="7"
map[29]="43"
map[30]="11"
map[114]="40"
map[32]="33"
map[33]="12"
map[34]="1862"
map[35]="45"
map[36]="H-34"

for i in 28 29 30 114 32 33 34 35 36 
do
    j=${map[$i]}

    conda activate dfastenv

    dfast -g ~/mycoplasma/genomes/${i}/assembly.fasta \
    --organism "Mycoplasma hominis" \
    --strain "H-34" \
    --locus_tag_prefix MHO_${j} \
    --minimum_length 200 \
    --references ~/mycoplasma/references/GCF_000085865.1_ASM8586v1_genomic.faa \
    --aligner blastp \
    --gcode 4 \
    --force \
    --out ~/mycoplasma/genomes/${i}/dfast_${i}

    conda activate ngsenv
    telegram-send dfast_${j}

done

wc -l ~/mycoplasma/genomes/*/dfast_*/proteins*

cat ~/mycoplasma/genomes/*/dfast_36*/protein.faa > all_proteins.faa
cat ~/mycoplasma/genomes/*/dfast_!(*36*)/protein.faa >> all_proteins.faa

srun --partition=RT --ntasks=1 --cpus-per-task=4 cd-hit -i all_proteins.faa -o clustered_proteins.faa -c 0.99 -n 5

# 1. Отфильтровываем кластеры
awk 'BEGIN{split("40|7|43|11|33|12|1862|45|H-34", ids, "|"); total_ids = length(ids);} />Cluster/{found_all = 1; for(i in seen) {if(seen[i] == 0) {found_all = 0; break;}} if(found_all && NR != 1) print s; s = $0; delete seen; for(i=1; i<=total_ids; i++) {seen[ids[i]] = 0;} next;} {for(id in ids) {if($0 ~ ids[id]) {seen[ids[id]] = 1;}} s = s "\n" $0;} END{found_all = 1; for(i in seen) {if(seen[i] == 0) {found_all = 0; break;}} if(found_all) print s;}' clustered_proteins.faa.clstr > filtered_clusters.faa.clstr
# 2. Извлекаем идентификаторы белков из отфильтрованных кластеров
grep '^0' filtered_clusters.faa.clstr | awk -F',' '{print $2}' | awk -F'|' '{print $1}' > selected_proteins.txt
# 3. Оставляем CDS
grep "CDS" genome_H-34.gff | awk -F';' '{for(i=1;i<=NF;i++) if ($i ~ /ID=/) print $i}' | cut -d'=' -f2 | grep -Fwf - selected_proteins.txt > filtered_proteins.txt
sed 's/^ >//' filtered_proteins.txt > corrected_filtered_proteins.txt

# 4. Оставляем только те белки, которые соответствуют выбранным идентификаторам
grep -A 1 -Fwf selected_proteins.txt protein_H-34.faa > core_proteins.faa