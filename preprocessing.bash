merge_fastq.pl Lk1_S1_R1_001.fastq Lk1_S1_R2_001.fastq Lk1_R1R2.fastq
#  use perl script written by 'simonandrews' in SEQanswers forum
# http://seqanswers.com/forums/showthread.php?t=10591
# extract UMIs
umi_tools extract --extract-method=string --bc-pattern=NNNNNN -L extract.log -I Lk1s.R1R2.fastq -S Lk1s.R1R2.exUMI.fastq
# split, trim, mask, align, annotate, dedup, count.
cat Lk1_R1R2_exUMI.fastq | fastx_barcode_splitter.pl --bcfile bc_index_file --prefix Lk1.R1R2.exUMI_ --bol --mismatches 1
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_007 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_007_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_011 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_011_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_012 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_012_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_013 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_013_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_014 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_014_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_015 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_015_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_016 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_016_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_017 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_017_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_018 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_018_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_019 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_019_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_020 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_020_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_021 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_021_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_022 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_022_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_023 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_023_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_024 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_024_masked
fastx_trimmer -Q 33 -f 20 -i Lk1.R1R2.exUMI_025 | fastq_masker -Q 33 -zv -o Lk1.R2.exUMI_025_masked
for i in Lk1.R2.exUMI*masked; do mkdir ${i/_masked/}; mv ${i} ./${i/_masked/}/${i}
cd ${i/_masked/}; STAR --genomeDir /home/maxwells/cel.seq/ref/dm6/ --sjdbGTFfile /home/maxwells/cel.seq/ref/dm6/dm6.refGene.CEL-Seq.gtf --readFilesIn ${i} --runThreadN 16 --readFilesCommand zcat --outReadsUnmapped Fastx
htseq-count -o ${i/_masked/}\.ann.sam -f sam Aligned.out.sam /home/maxwells/cel.seq/ref/dm6/dm6.refGene.CEL-Seq.gtf > ${i/_masked/}\.prededup.counts 
samtools view -H Aligned.out.sam > header
cat header ${i/_masked/}\.ann.sam > ${i/_masked/}\.ann.h.sam
rm ${i/_masked/}\.ann.sam
rm Aligned.out.sam
samtools view -bSt /home/maxwells/cel.seq/ref/dm6/dm6_ERCC_GAL4_GFP.fa.fai ${i/_masked/}\.ann.h.sam | samtools sort -m 4G -@ 10 -o ${i/_masked/}\.ann.bam
samtools index -@ 8 ${i/_masked/}\.ann.bam
umi_tools dedup --per-gene --gene-tag=XF --method=directional --log2stderr --output-stats=${i/_masked/}\.ann.dir -I ${i/_masked/}\.ann.bam -S ${i/_masked/}\.ann.dir.dedup.bam
htseq-count -f bam ${i/_masked/}\.ann.dir.dedup.bam /home/maxwells/cel.seq/ref/dm6/dm6.refGene.CEL-Seq.gtf > ${i/_masked/}\.dir.dedup.counts
mv ${i/_masked/}\.dir.dedup.counts ../
cd ..
done
