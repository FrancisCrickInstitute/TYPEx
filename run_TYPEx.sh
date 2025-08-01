#!/bin/bash
#LOAD MODULES
ml purge
ml Nextflow/22.04.0
ml Singularity/3.6.4
# export cache directory for singularity
export NXF_SINGULARITY_CACHEDIR='.cache'

release="test"
# Run TYPEx in standalone mode using the testdata data/cell_objects.tracerx.txt
nextflow run main.nf \
   -c $PWD/test.config \
   --input_dir $PWD/results/$release/ \
   --input_table $PWD/data/cell_objects.tracerx.txt \
   --sample_file $PWD/data/sample_file.tracerx.txt \
   --release $release \
   --output_dir "$PWD/results/$release/" \
   --params_config "$PWD/data/typing_params.json" \
   --annotation_config "$PWD/data/cell_type_annotation.testdata.json" \
   --overlay_config_file "$PWD/conf/tissue_segmentation.json" \
   --tissue_seg_model "$PWD/models/tumour_stroma_classifier.ilp" \
   --color_config $PWD/conf/celltype_colors.json \
   --deep_imcyto false --cellprofiler false \
   --exclude_cell_lineage "Vim+ cells" \
   -profile singularity \
   -resume
