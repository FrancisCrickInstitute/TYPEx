ml Nextflow/21.04.3
ml Singularity/3.6.4

# Set Nextflow singularity cache directory
# This folder is a symbolik link to a folder outside of the home directory, which has a limited quota
export NXF_SINGULARITY_CACHEDIR=$HOME/labwd/tools/singularity/cache
[[ ! -d $NXF_SINGULARITY_CACHEDIR ]] && mkdir $NXF_SINGULARITY_CACHEDIR

inputDIR="$PWD/data"
outDIR="$PWD/output"

study='schuerch'

nextflow run main.nf -profile singularity \
			-c nextflow.config \
			--inDir $inputDIR \
			--regFile $inputDIR/metadata.$study.txt \
			--inputTable cell_objects.schuerz.txt \
			-params-file conf/tissue_segmentation.json \
			--major_markers 'major_codex' \
			--stratify_label $study \
			--study $study \
			--outDir $outDIR \
			--imcyto false \
			--subsample false \
			-resume \
			-with-trace \
			-with-dag flowchart.png 

# -resume -with-report report.html \

