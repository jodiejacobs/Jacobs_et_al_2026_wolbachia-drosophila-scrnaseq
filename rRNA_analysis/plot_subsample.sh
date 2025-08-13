#!/bin/bash
#SBATCH --job-name=bam_sort
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=bam_sort_%j.out
#SBATCH --error=bam_sort_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jomojaco@ucsc.edu

# sbatch  plot_subsample.sh $OUTDIR

SCRIPT="/private/groups/russelllab/jodie/scRNAseq/wolbachia_titer/rRNA_analysis/scripts/plot_coverage.py"

OUTDIR=$1
SUB="0.0001 0.001 0.01 0.05 0.1 0.25 0.5"

for S in $SUB; do  
    echo "Plotting subsample $S"

    ls $OUTDIR/*coverage_${S}.tsv > $OUTDIR/coverage_files_${S}.txt
    python $SCRIPT $OUTDIR/coverage_files_${S}.txt --output-dir $OUTDIR/rRNA_plots_${S}

    echo "Completed plotting $S"
done



