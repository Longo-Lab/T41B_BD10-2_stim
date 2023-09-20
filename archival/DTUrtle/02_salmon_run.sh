#!/usr/bin/env bash

#SBATCH --job-name=DTU_pre2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=batch
#SBATCH --account=flongo
#SBATCH --time=24:00:00
#SBATCH --mem=128G


# load modules
ml anaconda htslib parallel

# conda environment (needed for salmon v 1.4.0)
source $(conda info --base)/etc/profile.d/conda.sh
conda activate salmon || { echo "salmon Conda environment not activated"; exit; }

# log the output of the run
exec 1>> $(date -d "today" +"%Y%m%d%H%M")_command.log 2>&1
set -ex

cd salmon
# getting reference sequences
SPEC="Mus_musculus"
BUILD="GRCm39"
VERSION="105"
URL="ftp://ftp.ensembl.org/pub/release-$VERSION"

# C wget -N $URL/fasta/${SPEC,}/dna/$SPEC.$BUILD.dna.toplevel.fa.gz
# C wget -N $URL/gtf/${SPEC,}/$SPEC.$BUILD.$VERSION.gtf.gz
# C wget -N $URL/fasta/${SPEC,}/cdna/$SPEC.$BUILD.cdna.all.fa.gz
# C wget -N $URL/fasta/${SPEC,}/ncrna/$SPEC.$BUILD.ncrna.fa.gz
# C # wget -N $URL/variation/vcf/${SPEC,}/${SPEC,}.vcf.gz
# C 
# C # naming variables
# C FASTA="$SPEC.$BUILD.dna.toplevel.fa.gz"
# C GTF="$SPEC.$BUILD.$VERSION.gtf.gz"
# C CDNA="$SPEC.$BUILD.cdna.all.fa.gz"
# C NCRNA="$SPEC.$BUILD.ncrna.fa.gz"
# C # VCF="${SPEC,}.vcf.gz"
# C # DBSNP="${SPEC,}.dbSNP140.vcf.gz"
# C 
# C # Building Salmon index
# C # make decoys list from genome
# C grep "^>" <(pigz -dc "$FASTA") | cut -d " " -f 1 | sed 's/>//g' > decoys.txt
# C # combine transcripts and decoys
# C cat $CDNA $NCRNA $FASTA > $BUILD.cdna.ncrna.decoy.fa.gz
# C rm $SPEC.$BUILD.cdna.all.fa.gz $SPEC.$BUILD.ncrna.fa.gz
# C # index
# C salmon index \
# C   -t $BUILD.cdna.ncrna.decoy.fa.gz \
# C   -d decoys.txt \
# C   -i $BUILD.quasi_index_23 \
# C   -k 23 \
# C   -p 16

cd ../fastq
parallel -j 2 salmon quant \
  -l A \
  -i ../salmon/$BUILD.quasi_index_23/ \
  -1 BD10_2_{}_1_trimmed.fastq.gz \
  -2 BD10_2_{}_2_trimmed.fastq.gz \
  --seqBias \
  --gcBias \
  --validateMappings \
  --rangeFactorizationBins 4 \
  --threads 8 \
  -o ../salmon/BD10_2_{}_trimmed.fq.gz/ \
  ::: 10 11 1 12 17 18 19 20 2 25 26 27 28 3 4 9

squeue -j $SLURM_JOBID --Format=TimeUsed
