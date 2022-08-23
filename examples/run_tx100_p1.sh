ml Nextflow/21.04.3
ml Singularity/3.6.4

export NXF_SINGULARITY_CACHEDIR=$HOME/labwd/tools/singularity/cache
[[ ! -d $NXF_SINGULARITY_CACHEDIR ]] && mkdir $NXF_SINGULARITY_CACHEDIR

panel='p1'
study='tracerx'
cohort='tx100'
nextflow run main.nf -profile singularity -resume -with-report report.html -with-trace -with-dag flowchart.png \
		-c nextflow.config --panel $panel --inputDir "/camp/project/proj-tracerx-lung/tctProjects/rubicon/${study}/${cohort}/imc/outputs/nextflow/"

# ingularity run  docker://mihangelova/typing:conda  perl -MFile::Find::Rule -e 'print "test\n"'
