# Npal_vir
## Sequencing data
ONT sequencing was performed for XX samples aiming to detect viruses and viroids.
## Processing
Options for next commands, change depending on local conditions:
```
threads=10
date=$(date +"%Y_%m_%d")
raw_files_pod5="path/to/file"
```
### Step 1: Basecalling with dorado
```
dorado basecaller sup $raw_files_pod5 --device cuda:0 --batchsize 64 > ${date}_${raw_files_pod5%pod5}bam 2> ${date}_out/${date}_dorado.log
```
Sort bam file and convert to fastq
```
samtools sort -@ $threads -n ${date}_${raw_files_pod5%pod5}bam -o ${date}_${raw_files_pod5%.pod5}_sorted.bam
```
Convert bam files to fastq
```
bedtools bamtofastq -i ${date}_${raw_files_pod5%.pod5}_sorted.bam -fq ${date}_${raw_files_pod5%.pod5}_sorted.fastq
```
### Step 2: Quality control
Do the following for each sample. Can be done iterating through a sample list, for example.
```
ls *.fastq | sed 's/.fastq//' > sample_list.txt
while read sample; do some_command ${sample}.fastq > ${sample}_processed.fastq; done < sample_list.txt
```
Get statistics before anything:
```
seqkit stats --all --tabular *.fastq > ${date}_stats_all.txt
```
Remove duplicates (by name first, reads duplicated for some reason), and then by sequence (exact clones).
```
seqkit rmdup --by-name $sample.fastq > ${sample}_rmdup.fastq
seqkit rmdup --by-seq ${sample}_rmdup.fastq -d ${sample}_duplicate_seqs.fastq -D ${sample}_duplicate_counts.txt -o ${sample}_rmdup2.fastq
```
Get statistics after deduplication:
```
seqkit stats --all --tabular *_rmdup2.fastq > ${date}_stats_rmdup.txt
```
Remove very short (length < 100) and poor-quality (q < 7) reads.
```
seqkit seq -m 100 -Q7 ${sample}_rmdup2.fastq > ${sample}_rmdup2_l100_q7.fastq
```
Get statistics after quality control:
```
seqkit stats --all --tabular *_rmdup2_l100_q7.fastq > ${date}_stats_rmdup_l100_q7.txt
```
Change lecture names for easier downstream interpretation.
```
seqkit fx2tab ${sample}_rmdup2_l100_q7.fastq | awk -v sample}="$sample" '{print sample "_" NR "\t" $0}' | cut -f1,3,4 | seqkit tab2fx > ${sample}_rmdup2_l100_q7_names.fastq
```
### Step 3: Host decontamination
Pool all reads maintaining source information to perform the decontamination step once and to increase information for clustering in Step 4.
```
cat *_rmdup2_l100_q7_names.fastq >> Pooled_rmdup2_l100_q7_names.fastq
```
Map to host (*Neltuma pallida*) available sequences (nuclear genome and organellar genomes).
```
minimap2 -ax map-ont ppa_v2.asm.fasta Pooled_rmdup2_l100_q7_names.fastq > Pooled_rmdup2_l100_q7_names_vsNG.sam
minimap2 -ax map-ont CP_Npallida.fasta Pooled_rmdup2_l100_q7_names.fastq > Pooled_rmdup2_l100_q7_names_vsCP.sam
minimap2 -ax map-ont MT_Npallida.fasta Pooled_rmdup2_l100_q7_names.fastq > Pooled_rmdup2_l100_q7_names_vsMT.sam
```
Get unmapped reads for each genome.
```
samtools -fastq -n -f 4 Pooled_rmdup2_l100_q7_names_vsNG.sam > Pooled_rmdup2_l100_q7_names_vsNG_f4.fastq
samtools -fastq -n -f 4 Pooled_rmdup2_l100_q7_names_vsCP.sam > Pooled_rmdup2_l100_q7_names_vsCP_f4.fastq
samtools -fastq -n -f 4 Pooled_rmdup2_l100_q7_names_vsMT.sam > Pooled_rmdup2_l100_q7_names_vsMT_f4.fastq
```
Select reads not mapping to any of the host sequences.
```
cat *_f4.fastq | seqkit fx2tab | cut -f1 | sort | uniq -c | sort -rn > Pooled_rmdup2_l100_q7_names_vsCPMTNG_f4.counts # Number of genomes to which each read DID NOT aligned to
grep " 3 " Pooled_rmdup2_l100_q7_names_vsCPMTNG_f4.counts | sed 's/ /\t/' | cut -f2 | seqkit grep --by-name -f - Pooled_rmdup2_l100_q7_names.fastq > Pooled_rmdup2_l100_q7_names_noCPMTNG.fastq
```
Prepare sequences for a blastn search. Create the database.
```
seqkit fq2fa Pooled_rmdup2_l100_q7_names_noCPMTNG.fastq > Pooled_rmdup2_l100_q7_names_noCPMTNG.fasta
makeblastdb -in RefSeq_TargetedLoci_Pcin_Nalb.fasta -dbtype nucl -out RefSeq_TargetedLoci_Pcin_Nalb.fasta
```
Search against a rRNA database using blastn. Used RefSeq Targeted Loci Project (2025-08-28) including Bacteria and Archaea 16S and 23S, and Fungi 18S, 28S and ITS. Also, added LSU and SSU and organellar rRNA from *P. cineraria* and *N. alba*.  
Notice the coverage (60%) and identity (60%) thresholds for this step and adjust if necessary.
```
blastn -query Pooled_rmdup2_l100_q7_names_noCPMTNG.fasta -db RefSeq_TargetedLoci_Pcin_Nalb.fasta -out Pooled_rmdup2_l100_q7_names_noCPMTNG.blastn -qcov_hsp_perc 60 -perc_identity 60 -outfmt 6 -num_threads $threads -max_hsps 2 -num_alignments 2
```
Get reads that did not match any rRNA sequence in the database.
```
cut -f1 Pooled_rmdup2_l100_q7_names_noCPMTNG.blastn | sort | uniq | seqkit grep --by-name -v -f Pooled_rmdup2_l100_q7_names_noCPMTNG.fastq > Pooled_rmdup2_l100_q7_names_noCPMTNG_noblast.fastq
```
### Step 4: Read clustering
Use a greedy, quality value-based algorithm for de novo clustering (10.1089/cmb.2019.0299). Use initial data in case the decontamination removed any viral read that may be similar to regions in the host genome or rRNA sequences.
```
isONclust --ont --fastq Pooled_rmdup2_l100_q7_names.fastq --outfolder isONclust_clusters/ --t 20
```
Also cluster the decontaminated dataset to get possible viral representatives.
```
isONclust --ont --fastq Pooled_rmdup2_l100_q7_names_noCPMTNG_noblast.fastq --outfolder Pooled_rmdup2_l100_q7_names_noCPMTNG_noblast_isONclust/ --t $threads
```
In both cases write fastq results for later analyses.
```
isONclust write_fastq --clusters isONclust_clusters/final_clusters.tsv --fastq Pooled_rmdup2_l100_q7_names.fastq --outfolder isONclust_seqs/ --N 1
isONclust write_fastq --clusters Pooled_rmdup2_l100_q7_names_noCPMTNG_noblast_isONclust/final_clusters.tsv --fastq Pooled_rmdup2_l100_q7_names_noCPMTNG_noblast.fastq --outfolder Pooled_rmdup2_l100_q7_names_noCPMTNG_noblast_isONclust_seqs/ --N 1
```
Get the cluster representatives of the clean data for searching against virus databases.
```
sed 's/\t/_/' Pooled_rmdup2_l100_q7_names_noCPMTNG_noblast_isONclust/final_cluster_origins.tsv | seqkit tab2fx > Pooled_rmdup2_l100_q7_names_noCPMTNG_noblast_isONclust_origins.fasta
```
Notice that the read name now contains the cluster number for the clean dataset, underscore, sample name, underscore and read number.
### Step 5: Read classification
For classification, the RefSeq Virus nucleotide release (2025-08) was used. A threshold of 50% coverage and identity was used.
```
makeblastdb -in viral.1.1.genomic.fna -dbtype nucl -out viral.1.1.genomic.fna
blastn -query Pooled_rmdup2_l100_q7_names_noCPMTNG_noblast_isONclust_origins.fasta -db viral.1.1.genomic.fna -out Pooled_rmdup2_l100_q7_names_noCPMTNG_noblast_isONclust_origins_viral.1.1.genomimc.fna.blastn -qcov_hsp_perc 50 -perc_identity 50 -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -num_threads $threads -max_hsps 5 -num_alignments 5
```
For representative reads matching viruses, get read counts for each sample based on the data before decontamination.
```
rep_read_name="read1" # Remove the cluster number for the clean dataset and underscore
unclean_cluster_n=$(grep -P "\t$rep_read_name$" isONclust_clusters/final_clusters.tsv | cut -f1)
grep -P "^$unclean_cluster_n\t isONclust_clusters/final_clusters.tsv | cut -f2 | cut -d'_' -f1 | sort | uniq -c | sort -rn # This removes the read number and leaves the sample name
```
Get the representative sequence for further analyses. For example, blastn against all GenBank sequences to discard obvious false positives, or further validate the results.
```
seqkit grep --by-name -p $rep_read_name Pooled_rmdup2_l100_q7_names.fastq
```
Additional databases were searched against (all downloaded on 2025/08/30):  
|DB name|type|num_seqs|sum_len|
|----|-----|-----|-----|
|viral.1.1.genomic.fna|Nucl|18,760|557,838,499|
|RefSeq Viral Protein non-human|Prot|638,295|149,664,160|
|RefSeq+GenBank Viral Protein non-human|Prot|7,898,515|2,131,752,193|
|RefSeq Viral Nucleotide non-human|Nucl|14,757|498,225,851|
|RefSeq+GenBank Viral Nucleotide non-human|Nucl|1,444,870|6,541,546,647|
|Viroids (viroids.org) All|Nucl|9,691|3,649,347|  

Draft: rustic code to get counts and representative sequence for all clusters with matches. For later embelishment:
```
while read seq_name; do vir_name=$(echo $seq_name | cut -d'_' -f2,3); new_cluster=$(grep -P "\t$vir_name$" ../isONclust_clusters/final_clusters.tsv | cut -f1); origin_seq=$(grep -P "^$new_cluster\t" ../isONclust_clusters/final_cluster_origins.tsv | cut -f3); printf "$seq_name\t$new_cluster\t"; while read sample; do count=$(grep -c -P "^$new_cluster\t$sample" ../isONclust_clusters/final_clusters.tsv); printf "$count\t"; done < ../sample_list.txt; printf "$origin_seq\n"; done < Pooled_rmdup2_l100_q7_names_noCPMTNG_noblast_isONclust_origins_RefSeq_Nucleotide.names | tee Pooled_rmdup2_l100_q7_names_noCPMTNG_noblast_isONclust_origins_RefSeq_Nucleotide.counts
```
Since there are many matches to make it unfeasible to look individually, all were searched locally against core-nt NCBI DB and visualized in MEGAN. Only representative sequences assigned to Virus or not assigned were considered.
## Software
BLAST v2.16.0+ - 10.1016/S0022-2836(05)80360-2  
minimap2 v2.24-r1122 - 10.1093/bioinformatics/bty191  
samtools v1.6 - 10.1093/bioinformatics/btp352  
isONclust v0.0.6.1 - 10.1089/cmb.2019.0299  
seqkit v2.8.2 - 10.1371/journal.pone.0163962  
dorado vXXXXX  
bedtools vXXXX
