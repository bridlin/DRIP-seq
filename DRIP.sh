#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type END
#SBATCH --mail-user b-barckmann@chu-montpellier.fr
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 6
#SBATCH --mem 64GB


module load cutadapt/4.0
module load fastqc/0.11.9
module load bowtie2/2.4.1
module load samtools/1.13
module load trim-galore/0.6.5
module load deeptools/3.5.0

input_list=("ERR2037035" \
	"ERR2037036" \
	"ERR2071902" \
	"ERR2071903" \
	"ERR2071906" \
	"ERR2071907")
input_dir=Briggs_2018-2019_wt-DRIP_fastq/
output_dir=Drip-seq_analysis
genome=genome/TriTrypDB-55_TbruceiLister427_2018_Genome/index_bowtie2/TriTrypDB-55_TbruceiLister427_2018_Genome
# genome=genome/TriTrypDB-62_TbruceiTREU927_Genome/index_bowtie2/TriTrypDB-62_TbruceiTREU927_Genome
prefix=aln_Tb427_2018_v55
# prefix=aln_Tb927_v62

# mkdir $output_dir
# for x in "${input_list[@]}"; do
# trim_galore --paired --output_dir $input_dir/ $input_dir/$x\_1.fastq.gz  $input_dir/$x\_2.fastq.gz && 
# # Briggs 2018 used the -very-sensitive option of bowtie2
# bowtie2 --threads 6 -x $genome -1 $input_dir/$x\_1_val_1.fq.gz -2 $input_dir/$x\_2_val_2.fq.gz   -S $output_dir/$x\_$prefix\.sam &&
# samtools view  -b -q 20 $output_dir/$x\_$prefix\.sam > $output_dir/$x\_$prefix\_filtered.bam &&
# samtools sort $output_dir/$x\_$prefix\_filtered.bam -o $output_dir/$x\_$prefix\_filtered_sorted.bam &&
# rm -f  $output_dir/$x\_$prefix\.sam &&
# rm -f  $output_dir/$x\_$prefix\_filtered.bam  ; done

input=("ERR2037035" \
	"ERR2071902" \
	"ERR2071906")
sample=("ERR2037036" \
	"ERR2071903" \
	"ERR2071907")

for y in 0 1 2 ; do
	bamCompare  --scaleFactorsMethod SES --binSize 50 --outFileFormat bedgraph -b1 $output_dir/${sample[$y]}\_$prefix\_filtered_sorted.bam -b2 $output_dir/${input[$y]}\_$prefix\_filtered_sorted.bam --outFileFormat bedgraph --outFileName $output_dir/${sample[$y]}.bdg ;done


