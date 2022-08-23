

process csm_submit {

	tag "CSM"
	label 'xs'
	executor = 'local'

//	container='mihangelova/typeconda:latest'

	//SCRIPTDIR=/camp/home/angelom/lung/scripts/csm
	///export PATH=/camp/home/angelom/labwd/tools/miniconda2/bin:$PATH
	//export PYTHONPATH=/camp/home/angelom/labwd/tools/miniconda2/bin/:$PYTHONPATH
	//sbatch -c 1 --ntasks=$NUMTHREADS --time=3-00:00:00  --mail-type=END --mail-user=mihaela.angelova@crick.ac.uk --mem-per-cpu=50G --wrap="python $SCRIPTDIR/cell_specific_masks.py $panel $run $cohort $study"

	input:
		tuple val(imageID), path(tifFile), path(inFile)
		val(csmDir)

	output:
		val (csmDir), emit: csmDir

	script:
	"""
		echo ${inFile} ${csmDir}
		cell_specific_masks_per_image.py ${inFile} ${tifFile} ${imageID} ${csmDir} 
		  # ${params.inDir}/${params.panel}/${params.run} ${params.outDir}/major/csm ${params.panel} ${params.run}
	"""
} 


process csm_export {
	
	tag "CSM export"
	label 'mem_medium'

	script:
	// sbatch --ntasks=1 --time=8:00:00  --mail-type=END \
//			--mail-user=mihaela.angelova@crick.ac.uk --mem-per-cpu=100G \
//			--wrap="Rscript $SCRIPTDIR/csm_exporter.R --run $run --panel $panel --cohort $cohort --study $study "

	input:
		val(csmDir)
		val files
	output:
		val params.outDir

	"""
		export BASE_DIR=$baseDir
		csm_exporter.R --inDir ${csmDir} \
			--run ${params.run} --panel ${params.panel} --cohort ${params.cohort} \
			--study ${params.study} --celltypeReviewFile ${params.celltypeReviewFile} \
			--markers ${params.major_markers} --nndist 5 \
			--tissAreaDir "${params.outDir}/tissue_seg/" 
	"""
}
