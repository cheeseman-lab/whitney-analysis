#!/bin/bash
# Run SBS processing for plate 1 using config_plate1.yml

mkdir -p slurm/slurm_output/main
start_time_formatted=$(date +%Y%m%d_%H%M%S)
log_file="slurm/slurm_output/main/sbs-plate1-${start_time_formatted}.log"
exec > >(tee -a "$log_file") 2>&1

start_time=$(date +%s)
echo "===== STARTING SBS PROCESSING FOR PLATE 1 ====="
echo "Using config: config/config_plate1.yml"

snakemake --use-conda --cores 100 \
    --snakefile "../brieflow/workflow/Snakefile" \
    --configfile "config/config_plate1.yml" \
    --rerun-triggers mtime \
    --rerun-incomplete \
    --keep-going \
    --groups align_sbs=sbs_tile_group \
            apply_ic_field_sbs=sbs_tile_group \
            segment_sbs=sbs_tile_group \
            extract_sbs_info=sbs_tile_group \
    --until all_sbs \
    --config plate_filter=1

if [ $? -ne 0 ]; then
    echo "ERROR: Processing of plate 1 failed."
fi

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "===== PLATE 1 COMPLETED ====="
echo "Runtime: $((duration / 3600))h $(((duration % 3600) / 60))m $((duration % 60))s"
