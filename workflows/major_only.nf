
///////////////////////////////////////////////
/* --        INCLUDE SUBWORKFLOWS        -- */
//////////////////////////////////////////////

include { format_input; collate_features } from '../modules/format.nf'
include { csm_submit; csm_export } from '../modules/csm.nf'
include { call_major; call_subtypes; call_major_ref; review_major_types; call_subtypes_nostrat } from '../modules/typing.nf'
include { preprocess_panck; create_composites; run_classifier; process_probabilities; overlay; export_tissue_seg } from '../modules/tissue_seg.nf'
include { call_method; call_subsampled; match_clusters; plot_subsampled } from '../modules/robustness.nf'
include { exporter; exporter_major; plot_dr } from '../modules/release.nf'
include { compare_dp; compare_subsampled; compare_tcra; compare_panels; compare_flow } from '../modules/summary.nf'


features = Channel.of( ['MeanIntensity'], ['LocationCenter'] , ['AreaShape'])
subsample_methods = Channel.of(['flowSOM'], ['Rphenograph'], ['FastPG'],  ['cellassign']  )
subsample_markers = Channel.of( ['mcsa'] )
iterations = Channel.of( [1], [2], [3] )
subtype_methods_nostrat = Channel.of(['Rphenograph_FALSE'], ['FastPG_FALSE'], ['flowSOM_FALSE'] )
subtype_methods = Channel.of(['Rphenograph'], ['FastPG'], ['flowSOM'] )


workflow MAJOR_ONLY {

	if(params.imcyto) {

		// what if they are not named run
		cellFilePattern="${params.inDir}/${params.panel}/${params.run}/run*/results/segmentation/**/Cells.csv"
		samplePattern=cellFilePattern.replaceAll('\\*', '.*').replaceAll('.\\*.\\*', '(.*)').replaceAll('\\/+', '.')
	
		cellFiles=Channel
			.fromPath("${cellFilePattern}", relative:false)
			.map{ file -> tuple(file.toString().replaceAll(samplePattern, '$1').replaceAll('\\/', '-'), file) }

		format_input(cellFiles, "${params.outDir}/nfData")
		collate_features("${params.outDir}/nfData/", "${params.outDir}/features", features, 
			format_input.out.collect{ it.baseName }.unique())
		
		rawMasksPattern=cellFilePattern="${params.inDir}/${params.panel}/${params.run}/run*/results/segmentation/**/raw*.tiff"
		rawSamplePattern=rawMasksPattern.replaceAll('\\*', '.*').replaceAll('.\\*.\\*', '(.*)').replaceAll('\\/+', '.')
		rawMasks=Channel
            .fromPath("${rawMasksPattern}", relative:false)
        	.map{ file -> tuple(file.toString().replaceAll(rawSamplePattern, '$1').replaceAll('\\/', '-'), file) }
		
		csm_submit(rawMasks.combine(cellFiles, by: 0), "${params.outDir}/major/${params.major_markers}/csm")
		csm_export("${params.outDir}/major/${params.major_markers}/csm", csm_submit.out.collect())
		
		// Requires format_input to have finalised
		call_major(params.major_method, "${params.outDir}/nfData", collate_features.out.collect{it.baseName}.unique())
		call_major_ref(params.major_method, "${params.outDir}/nfData", collate_features.out.collect{it.baseName}.unique())
	}

	 if(params.subsample) {
		call_method(subsample_methods, "${params.outDir}/nfData", collate_features.out.collect{it.baseName}.unique())
		exporter_major(subsample_methods, call_method.out.collect())

		pars=iterations.combine(subsample_methods).combine(subsample_markers)
		pars.view()
		call_subsampled(pars, "${params.outDir}/nfData", collate_features.out.collect{it.baseName}.unique())
		match_clusters(params.outDir, subsample_methods, call_subsampled.out.collect(), call_method.out.collect())
        // plot_subsampled(match_clusters.out.collect())
		plot_dr('umap', "${params.outDir}/nfData", call_method.out.collect())
    }

	preprocess_panck(params.inDir)
	create_composites(params.inDir, preprocess_panck.out)
	run_classifier(create_composites.out)
	process_probabilities(run_classifier.out)
	overlay(process_probabilities.out)
	export_tissue_seg(overlay.out)
	
	//	if(params.major_only == false) {	
	// do review - stratify
	review_major_types(call_major.out, call_major_ref.out, csm_export.out.collect(), export_tissue_seg.out)
	// allow to be optional
	
	call_subtypes(subtype_methods, "${params.outDir}/nfData", review_major_types.out, export_tissue_seg.out)
	call_subtypes_nostrat(subtype_methods, "${params.outDir}/nfData", review_major_types.out, export_tissue_seg.out)
	
	subtype_runs=subtype_methods.mix(subtype_methods_nostrat)
	subtype_runs.view()
	exporter(subsample_methods_nostrat, call_subtypes.out, export_tissue_seg.out, call_subtypes_nostrat.out)

	//	}
}
