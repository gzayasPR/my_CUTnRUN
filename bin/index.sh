source project_env.sh
ecoli_ref=${my_data}/references/E_Coli_K12_MG1655
cd ${ecoli_ref}
tar -xvzf Escherichia_coli_K_12_MG1655_NCBI_2001-10-15.tar.gz

mus_ref=${my_data}/references/Mus_musculus_Ensembl_GRCm38
cd ${mus_ref}
tar -xvzf Mus_musculus_Ensembl_GRCm38.tar.gz

# Define root directory
GENOME_DIR="${my_data}/references/Mus_musculus_GRCm39/Ensembl/111"

# Create the structure
mkdir -p ${GENOME_DIR}/Annotation/Genes
mkdir -p ${GENOME_DIR}/Sequence/{Bowtie2Index,Chromosomes,WholeGenomeFasta}
# Navigate to base directory
cd ${GENOME_DIR}

# Download genome FASTA
wget ftp://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz -P Sequence/WholeGenomeFasta/

# Download annotation GTF
wget ftp://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.gtf.gz -P Annotation/Genes/

# Decompress genome FASTA
gunzip Sequence/WholeGenomeFasta/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

gunzip Annotation/Genes/Mus_musculus.GRCm39.111.gtf.gz
awk '$3 == "transcript" {
  match($0, /gene_name "([^"]+)"/, gn);
  match($0, /transcript_id "([^"]+)"/, tid);
  gene_name = (gn[1] != "") ? gn[1] : tid[1];
  start = $4 - 1;
  end = $5;
  strand = $7;
  print $1"\t"start"\t"end"\t"gene_name"\t0\t"strand;
}' Annotation/Genes/Mus_musculus.GRCm39.111.gtf > Annotation/Genes/Mus_musculus.GRCm39.111_transcripts.bed


awk '$3 == "gene" {
  match($0, /gene_name "([^"]+)"/, gn);
  match($0, /gene_id "([^"]+)"/, gid);
  gene_name = (gn[1] != "") ? gn[1] : gid[1];
  start = $4 - 1;
  end = $5;
  strand = $7;
  print $1"\t"start"\t"end"\t"gene_name"\t0\t"strand;
}' Annotation/Genes/Mus_musculus.GRCm39.111.gtf > Annotation/Genes/Mus_musculus.GRCm39.111_genes.bed

# # Generate TSS BED file from GTF
# awk '$3 == "transcript" {
#   match($0, /gene_name "([^"]+)"/, gn);
#   match($0, /transcript_id "([^"]+)"/, tid);
#   gene_name = (gn[1] != "") ? gn[1] : tid[1];
#   if ($7 == "+") {
#     print $1"\t"($4 - 1)"\t"$4"\t"gene_name"\t0\t+"
#   } else if ($7 == "-") {
#     print $1"\t"($5 - 1)"\t"$5"\t"gene_name"\t0\t-"
#   }
# }' Annotation/Genes/Mus_musculus.GRCm39.111.gtf \
# > Annotation/Genes/Mus_musculus.GRCm39.111_TSS.bed

ml bowtie2
# Create index with Bowtie2
bowtie2-build Sequence/WholeGenomeFasta/Mus_musculus.GRCm39.dna.primary_assembly.fa \
              Sequence/Bowtie2Index/Mus_musculus.GRCm39
echo "Mouse genome (GRCm39) from Ensembl release 111" > README.txt
mv README.txt Annotation/

# Navigate to base directory
cd ${GENOME_DIR}
mkdir BLACKLIST
cd BLACKLIST


source project_env.sh
EC_FA="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/data/references/E_Coli_K12_MG1655/Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Sequence/Bowtie2Index/genome.fa"
MM_FA="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/data/references/Mus_musculus_GRCm39/Ensembl/111/Sequence/WholeGenomeFasta/Mus_musculus.GRCm39.dna.primary_assembly.fa"
# Modules (HiPerGator style)
module load bowtie2 samtools

mkdir -p ${my_data}/references/combined_ref/
cd ${my_data}/references/combined_ref/

# 1a) Prefix E. coli contig names so they never collide with mouse
# (keeps everything else intact)
EC_FA_TAGGED=Ecoli_K12.tagged.fa
awk 'BEGIN{OFS=""} /^>/ {sub(/^>/,">ecoli|"); gsub(/[ \t].*$/,""); print; next} {print}' \
  "$EC_FA" > "$EC_FA_TAGGED"

# 1b) Concatenate into one FASTA
cat "$MM_FA" "$EC_FA_TAGGED" > mouse_ecoli.concat.fa

# 1c) Build Bowtie2 index once
bowtie2-build mouse_ecoli.concat.fa mouse_ecoli
