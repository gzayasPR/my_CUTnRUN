# General settings
NAME="Run2"
account="mateescu"
cpus=8

source project_env.sh
proj_env=${my_bin}/project_env.sh

# Metadata
CSV="${my_data}/metadata/data.csv"

# References
ECOLI_INDEX="${my_data}/references/E_Coli_K12_MG1655/Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Sequence/Bowtie2Index/genome"
MOUSE_INDEX="${my_data}/references/Mus_musculus_GRCm39/Ensembl/111/Sequence/Bowtie2Index/Mus_musculus.GRCm39"
MOUSE_GTF="${my_data}/references/Mus_musculus_GRCm39/Ensembl/111/Annotation/Genes/Mus_musculus.GRCm39.111.gtf"
BLACKLIST="${my_data}/references/Mus_musculus_GRCm39/Ensembl/mm39.excluderanges/bed"

# ─── Step control flags (1/true/yes to run) ───────────
RUN_STEP1=0      # QC KMET
RUN_STEP2=0   # Alignment
RUN_STEP3=0      # Normalization
RUN_STEP4=0     # Peak calling
RUN_STEP5=1    # Annotation