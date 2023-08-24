mixcr analyze -f -Xmx30g milab-human-rna-tcr-umi-race \
    --tag-pattern "^GGGTCAGGGTTCTGGATAT(R1:*)\^(UMI:N{4}TN{4}TN{4})CTTG{:5}(R2:*)" \
    A/V350171692_Run204_L01_8_1.fq.gz \
    A/V350171692_Run204_L01_8_2.fq.gz \
    results/K_TRA 

mixcr analyze -f -Xmx30g milab-human-rna-tcr-umi-race \
    --tag-pattern "^AACACcTTtTTCAGGTCCT(R1:*)\^(UMI:N{4}TN{4}TN{4})CTTG{:5}(R2:*)" \
    B/V350171692_Run204_L02_8_1.fq.gz \
    B/V350171692_Run204_L02_8_2.fq.gz \
    results/K_TRB 


mixcr assemble --report results/K_TRA.assemble.report.txt --json-report results/K_TRA.assemble.report.json results/K_TRA.refined.vdjca results/K_TRA.clns

####################################################
#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate TCRenv
mixcr analyze generic-amplicon-with-umi \
    -f -Xmx30g \
    --species hsa \
    --rna \
    --rigid-left-alignment-boundary \
    --floating-right-alignment-boundary C \
    --tag-pattern "^AACACSTTKTTCAGGTCCT(R1:*)\^(UMI:N{4}TN{4}TN{4})CTTG{:5}(R2:*)" \
    B/V350171692_Run204_L02_8_1.fq.gz \
    B/V350171692_Run204_L02_8_2.fq.gz \
    results/K_TRB 

conda activate ngsenv
telegram-send TRB

conda activate TCRenv
mixcr analyze generic-amplicon-with-umi \
    -f -Xmx25g \
    --species hsa \
    --rna \
    --rigid-left-alignment-boundary \
    --floating-right-alignment-boundary C \
    --tag-pattern "^GGGTCAGGGTTCTGGATAT(R1:*)\^(UMI:N{4}TN{4}TN{4})CTTG{:5}(R2:*)" \
    A/V350171692_Run204_L01_8_1.fq.gz \
    A/V350171692_Run204_L01_8_2.fq.gz \
    results/K_TRA  

conda activate ngsenv
telegram-send TRA

mixcr refineTagsAndSort -Xmx16g --report results/K_TRA.refine.report.txt --json-report results/K_TRA.refine.report.json results/K_TRA.vdjca results/K_TRA.refined.vdjca

conda activate ngsenv
telegram-send --file alignQc.pdf



conda activate TCRenv
mixcr analyze generic-amplicon-with-umi \
    -f \
    --species hsa \
    --rna \
    --rigid-left-alignment-boundary \
    --floating-right-alignment-boundary C \
    --tag-pattern "^(R1:*)\^(UMI:N{14})CTTG{3:5}(R2:*)" \
    A/V350171692_Run204_L01_8_1.fq.gz \
    A/V350171692_Run204_L01_8_2.fq.gz \
    results/K_TRA  

conda activate ngsenv
telegram-send TRA

conda activate TCRenv
mixcr analyze generic-amplicon-with-umi \
    -f \
    --species hsa \
    --rna \
    --rigid-left-alignment-boundary \
    --floating-right-alignment-boundary C \
    --tag-pattern "^(R1:*)\^(UMI:N{14})CTTG{3:5}(R2:*)" \
    B/V350171692_Run204_L02_8_1.fq.gz \
    B/V350171692_Run204_L02_8_2.fq.gz \
    results/K_TRB 
    
conda activate ngsenv
telegram-send TRB


###/storage/data2/sikamov/TCR/results


####QC
fastqc B/*fq.gz -o qc/
telegram-send --file qc/*.html


sudo cp -r results/ /mnt/yandex.disk/Downloads/TCR_seq/


seqtk sample -s100 V350171692_Run204_L01_8_1.fq.gz 5279075 > V350171692_Run204_L01_8_1_sub.fq
seqtk sample -s100 V350171692_Run204_L01_8_2.fq.gz 5279075 > V350171692_Run204_L01_8_2_sub.fq

/Volumes/SATECHI/NGS/


conda activate TCRenv

mixcr analyze generic-amplicon-with-umi \
    -f -Xmx30g \
    --species hsa \
    --rna \
    --rigid-left-alignment-boundary \
    --floating-right-alignment-boundary C \
    --tag-pattern "^(R1:*)\^(UMI:N{14})CTTG{3:5}(R2:*)" \
    A/V350171692_Run204_L01_8_1_sub.fq \
    A/V350171692_Run204_L01_8_2_sub.fq \
    results/K_TRA_sub  

conda activate ngsenv
telegram-send K_TRA_sub
telegram-send --file results/K_TRA_sub*.txt




###################

srun --partition=RT --ntasks=1 --cpus-per-task=4 seqtk sample -s100 V350171692_Run204_L01_8_1.fq 5279075 > V350171692_Run204_L01_8_1_sub.fq

srun --partition=RT --ntasks=1 --cpus-per-task=4 seqtk sample -s100 V350171692_Run204_L01_8_2.fq 5279075 > V350171692_Run204_L01_8_2_sub.fq


#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=RT  # Замените YOUR_PARTITION на нужное значение
#SBATCH --job-name=TCR
#SBATCH --mail-user=sikamov.kv@phystech.edu  # Замените на вашу почту
#SBATCH --comment="TCR"

source ~/miniconda3/etc/profile.d/conda.sh 


conda activate TCRenv
mixcr analyze generic-amplicon-with-umi \
    -f \
    --species hsa \
    --rna \
    --rigid-left-alignment-boundary \
    --floating-right-alignment-boundary C \
    --tag-pattern "^(R1:*)\^(UMI:N{14})CTTG{3:5}(R2:*)" \
    B/V350171692_Run204_L02_8_1.fq.gz \
    B/V350171692_Run204_L02_8_2.fq.gz \
    results/K_TRB 
    
conda activate ngsenv
telegram-send TRB


conda activate TCRenv
mixcr analyze generic-amplicon-with-umi \
    -f \
    --species hsa \
    --rna \
    --rigid-left-alignment-boundary \
    --floating-right-alignment-boundary C \
    --tag-pattern "^(R1:*)\^(UMI:N{14})CTTG{3:5}(R2:*)" \
    A/V350171692_Run204_L01_8_1_sub.fq \
    A/V350171692_Run204_L01_8_2_sub.fq \
    results/K_TRA  

conda activate ngsenv
telegram-send TRA