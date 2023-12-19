source /home/d.omelchenko/mambaforge/bin/activate ngs_course

Reads=/home/m.logacheva/project/mitochondrial_disfunction/source_data
Genome=/home/m.logacheva/project/mitochondrial_disfunction/hg38_ref
Work=/home/m.logacheva/project/mitochondrial_disfunction_work

#Perform fastqc
mkdir fastqc_out

for file in $Reads*.fastq;
do
	fastqc -t 6 -o $Work/fastqc_out $file
done

#Perform trimming
mkdir trimmed_reads
mkdir trim_report

for file in $Reads/*.fastq;
do
	name=$(basename "$file")
	filename=${name%.*}
	fastp --in1 $file --out1 $Work/trimmed_reads/trimmed_${name} --cut_front --cut_right --thread 6 -h $Work/trim_report/${filename}_report.html
done

#Perform fastqc after trimming
mkdir trimmed_fastqc_out

for file in $Work/trimmed_reads/*.fastq;
do
	fastqc -t 6 -o $Work/trimmed_fastqc_results $file
done

#Index reference genome
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir genome_index --genomeFastaFiles $Genome/hg38.fa --sjdbGTFfile $Genome/hg38.ensGene.gtf --sjdbOverhang 49

#Mapping reads on reference genome
mkdir mapped

for file in $Work/trimmed_reads/*.fastq;
do
	name=$(basename "$file")
	filename=${name%.*}
	STAR --runThreadN 6 --genomeDir $Work/genome_index --readFilesIn $file --outFileNamePrefix $Work/mapped/${filename}.sam
done

#Count matrix
featureCounts -T 6 -t exon -g gene_id -a $Genome/hg38.ensGene.gtf -o reads/readcounts.txt $Work/mapped/*.sam

echo 'Script is done!'

