nextflow run main.nf -resume -profile singularity,test \
    --outdir results_make_targets \
    --blast_unassigned_sequences false \
    --make_targets true \
    --surveillance_columns refseq/surveillance_columns_research_make_targets.txt \
    --primer_mix iPM-01 \
    --study_type research \
    --all_refseq refseq/refseq_all.xlsx \