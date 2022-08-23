///////////////////////////////////////////////
/* --        INCLUDE SUBWORKFLOWS        -- */
//////////////////////////////////////////////
// TODO
// + process user-provided input without imcyto
// TODO: + PROCESS INPUT GENERAL OR IMCYTO
// TODO: + do not run nostrat if not required [set up params.stratify_by_confidence as boolean]
// TODO: * overlay based on tissue segmentation config json file
// TODO: * Verify input files/matrix are provided
// TODO: * TYPE as for subsampling methods [iteration as optional input]
// TODO: * remove call_method, call_subsampled and other from robustness
// TODO: * W or W/o imcyto
// TODO: * see how to pass on subsample
// TODO: * also exporter major in subsample not run
// TODO: * major only w/o the two-tiered process
// TODO: * Better place for cellObjFile?

// Lower priority
// TODO: LP: W or W/o raw tissue images 
// TODO: LP: HOW TO PARALELISE AND SPEED UP THE TISSUE SEG PROCESS; plus include input files
// TODO: LP: WHAT HAPPENS IF NO RAW IMAGES PROVDIDED? ONLY CHECK TISSSEG.OUT WHEN TISSEGDIR EXISTS ->
	 // for this tissue segmentation markers need to be provded from file, but now hardcoded
// TODO: LP: Should run even when TISSEG.out hasn't been finished in cases when the dir exists with user-provided masks
// TODO: LP: // plot_subsampled(match_clusters.out.collect())
// TODO: LP: make gpu server preferences for tsne
// TODO: LP: DO I need cellAssignReassign file as input
		
/////////////////////////////////////////////
// check that nostrat is not run on the same model;
// e.g. if not specified by the user or if a stratificiation model was not made
// i.e. manually set which marker to be excluded and check if present or on top of the tree?
		
/* -- LOAD MODULES and FUNCTIONS -- */
include { format_input; collate_features } from '../modules/format.nf'
include { csm_submit; csm_export } from '../modules/csm.nf'
include { TYPE as tier_one; TYPE as tier_two; TYPE as tier_one_ref; 
		  review_major_types; TYPE as tier_two_nostrat }  from '../modules/typing.nf'
include { preprocess_panck; create_composites; run_classifier; 
		  process_probs; ts_exporter; mask_overlay }  from '../modules/tissue_seg.nf'

include { TYPE as call_cluster } from '../modules/typing.nf'
include { call_subsampled; match_clusters; plot_subsampled }  from '../modules/robustness.nf'
include { exporter as subtypes_exporter; exporter as major_exporter; 
			plot_dr } from '../modules/release.nf'
include { get_cell_files; get_imcyto_raw_masks; parse_json_file;
			get_tissue_masks_config }  from '../lib/functions.nf'

/* -- INPUT TYPING PARAMETERS -- */
subtype_methods = 
	Channel.of(['FastPG']
			   //['Rphenograph'],
			   //['flowSOM']
			  )
subsample_methods =
	Channel.of(['flowSOM'],
			   ['Rphenograph'],
			   ['FastPG'],
			   ['cellassign']
			  )
features = 
	Channel.of(['MeanIntensity'],
			   ['LocationCenter'], 
			   ['AreaShape']
			  )
subsample_markers = Channel.of( [params.major_markers] )
iterations = Channel.of( [1], [2], [3] )

cellObjFile="${params.outDir}/features/LocationCenter/${params.panel}_LocationCenter_${params.run}.csv"

workflow TYPEx {
	
	/* -- INPUT PREPROCESSING -- */
	PREPROCESS()
	
	/* -- TISSUE SEGMENTATION -- */
	TISSEG()
	
	/* -- OVERLAY TISSUE ANNOTATIONS */
	overlayParams = get_tissue_masks_config('conf/tissue_segmentation.json')
	mask_overlay(
		overlayParams,
		cellObjFile,
		TISSEG.out,
		PREPROCESS.out.features
	)
	
	if(params.tiered) {
		println "Running tiered typing"
		TIERED(PREPROCESS.out)
	}

	if(params.sampled || params.clustered) {
		println "Calling subsampling"
		SUBSAMPLING(PREPROCESS.out)
	}
	
}

workflow PREPROCESS {

	main:
		cellFiles=get_cell_files(params.imcyto)
		// if(cellFiles.size() > 0) {
			cellFiles.combine(features).view()
			format_input(cellFiles.combine(features))
			collate_features(features, format_input.out.collect())
		// }
	emit:
		features=collate_features.out.collect()
}

workflow TISSEG {

	///////////////////////////////////////////////
	/* --        TISSUE SEGMENTATON        -- */
	//////////////////////////////////////////////
	
	main:
		preprocess_panck()
		create_composites(preprocess_panck.out)
		run_classifier(create_composites.out)
		process_probs(run_classifier.out)
		ts_exporter(process_probs.out)
	emit:
		ts_exporter.out
}

workflow STRATIFY {

	take: out
	main:
		cellFiles=get_cell_files(params.imcyto)
		rawMasks=get_imcyto_raw_masks()
		csm_submit(
			rawMasks.combine(cellFiles, by: 0), 
			"${params.outDir}/major/${params.major_markers}/csm"
		)
	
		csm_export(
			"${params.outDir}/major/${params.major_markers}/csm",
			csm_submit.out.collect())

		tier_one_ref(
			params.major_method,
			'major',
			params.major_markers,
			false,
			"${params.outDir}/nfData", 
			out.features
		)
	
		review_major_types(
			tier_one.out, 
			tier_one_ref.out,
			csm_export.out,
			TISSEG.out
		)		 		
	emit:
		review_major_types.out
}

workflow TIERED {

	take: out
	main:
		/* --        TIER 1        -- */
		tier_one(
			params.major_method,
			'major',
			params.major_markers,
			false,
			"${params.outDir}/nfData",
			out.features
		)
		if(params.imcyto && params.stratify_by_confidence) {
			// raw segmentation per celltype-specific markers are needed
			STRATIFY()
		}
	
		/* --        TIER 2        -- */	
		tier_two(
				subtype_methods,
				'subtypes',
				"${params.subtype_markers}",
				params.stratify_by_confidence,
				"${params.outDir}/nfData",
				STRATIFY.out && TISSEG.out
			)
		if(params.stratify_by_confidence) {
				subtypes_exporter(subtype_methods, 'subtypes', tier_two.out)
		} else {
			subtype_methods_nostrat = subtype_methods.map{ method -> [ "${method[0]}_FALSE" ]}
			subtypes_exporter(subtype_methods_nostrat, 'subtypes', tier_two.out)
		}
	emit:
		subtypes.explorer.out
}

workflow SUBSAMPLING {

	///////////////////////////////////////////////
	/* --        ROBUSTNESS ANALYSES        -- */
	//////////////////////////////////////////////
	take: out
	main:
		pars=iterations.combine(subsample_methods)
				.combine(subsample_markers)
		pars.view()
		call_cluster(
			subsample_methods, 
			"major", 
			params.major_markers,
			"false",
			"${params.outDir}/nfData",
			out
		)
	 	major_exporter(
			subsample_methods,
			'major',
			call_cluster.out.collect())
			
		call_subsampled(
            pars,
            "${params.outDir}/nfData",
			out
		)

		match_clusters(
			params.outDir, 
			subsample_methods, 
			call_subsampled.out.collect(), 
			call_cluster.out.collect()
		)

		plot_dr(
			"${params.dimred_method}", 
			"${params.outDir}/nfData", 
			call_cluster.out.collect()
		)
}
