///////////////////////////////////////////////
/* --        INCLUDE SUBWORKFLOWS        -- */
//////////////////////////////////////////////

		
/////////////////////////////////////////////
/* -- LOAD MODULES and FUNCTIONS -- */
include { format_input; collate_features } from '../modules/format.nf'
include { csm_submit; csm_export } from '../modules/csm.nf'
include { TYPE as tier_one; TYPE as tier_two; TYPE as tier_one_ref; 
		  TYPE as tier_two_nostrat; build_strata_model }  from '../modules/typing.nf'
include { create_composites; run_classifier; 
		  process_probs; ts_exporter; mask_overlay }  from '../modules/tissue_seg.nf'

include { TYPE as call_cluster } from '../modules/typing.nf'
include { call_subsampled; match_clusters; plot_subsampled }  from '../modules/robustness.nf'
include { exporter as subtypes_exporter; exporter as major_exporter; 
			plot_dr } from '../modules/release.nf'
include { get_cell_files; get_imcyto_raw_masks; parse_json_file;
		get_tissue_seg_markeres; get_tissue_masks_config }  from '../lib/functions.nf'
include { QC as QC_subtype; QC as QC_major } from '../workflows/qcw.nf'

		
/* -- INPUT TYPING PARAMETERS -- */
subsample_methods =
	Channel.of(
			 ['FastPG'],
			 ['Rphenograph'] 
			 // ['flowSOM'], ['cellassign']
	)
features = 
	Channel.of(['MeanIntensity'],
			   ['LocationCenter'],
			   ['AreaShape']
			  )
subsample_markers = Channel.of( [params.annotate_markers] )
iterations = Channel.of( [1], [2], [3] )

cellObjFile="${params.output_dir}/features/LocationCenter/${params.panel}_LocationCenter_${params.release}.csv"

workflow TYPEx {

	/* -- INPUT PREPROCESSING -- */
	PREPROCESS()
	
	/* -- TISSUE SEGMENTATION -- */
	TISSEG()
	
	/* -- OVERLAY TISSUE ANNOTATIONS */
	overlayParams = get_tissue_masks_config(params.overlay_config_file)
	mask_overlay(
		overlayParams,
		cellObjFile,
		TISSEG.out,
		PREPROCESS.out.features
	)
	
	if(params.tiered) {
		println "Running tiered typing"
		TIERED(PREPROCESS.out, mask_overlay.out.collect())
	}

	if(params.sampled || params.clustered) {
		println "Calling subsampling"
		SUBSAMPLING(PREPROCESS.out)
	}
	
}

workflow PREPROCESS {

	main:
		cellFiles=get_cell_files(params.deep_imcyto, params.cellprofiler)
		// if(cellFiles.size() > 0) {
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
	
		run_dirs = Channel.fromPath(params.image_dir, type: "dir", relative: false)
		tissegMarkers = get_tissue_seg_markeres(params.overlay_config_file)	
		create_composites(tissegMarkers, run_dirs)
		run_classifier(create_composites.out.collect())
		process_probs(run_classifier.out)
		ts_exporter(process_probs.out)
	emit:
		ts_exporter.out
}

workflow TIERED {

	take: out
		  mask_overlay
	main:
	
		/* --        TIER 1        -- */
		print("Tier 1")
		
		tier_one(
			params.major_method,
			'major',
			params.major_markers,
			params.major_markers,
			false,
			"${params.output_dir}/nfData",
			mask_overlay // waits for PREPROCESS, TISSEG and mask_overlay
		)
		if(params.stratify_by_confidence) {
			// raw segmentation per celltype-specific markers are needed
			println 'Stratifying and tier 2'
			STRATIFY(out, mask_overlay, tier_one.out)
	
			/* --        TIER 2        -- */	
			tier_two(
					params.subtype_method,
					'subtypes',
					"${params.subtype_markers}",
					"${params.major_markers}",
					params.stratify_by_confidence,
					"${params.output_dir}/nfData",
					STRATIFY.out
				)
		} else {
		
			println 'Tier 2 w/o stratification'
			/* --        TIER 2        -- */
            tier_two(
                    params.subtype_method,
                    'subtypes',
                    "${params.subtype_markers}",
                    "${params.major_markers}",
                    params.stratify_by_confidence,
                    "${params.output_dir}/nfData",
                    tier_one.out
                )

		}
		println "Exporter"
		if(params.stratify_by_confidence) {
			subtypes_exporter(params.subtype_method, 
				'subtypes',
				params.subtype_markers,
				params.major_markers,
				tier_two.out)
			QC_subtype(subtypes_exporter.out,
				"${params.subtype_markers}",
				'subtypes',
                params.subtype_method)
		} else {
			params.subtype_method_nostrat = params.subtype_method.map{ method -> [ "${method[0]}_FALSE" ]}
			subtypes_exporter(
				params.subtype_method_nostrat,
				'subtypes',
				params.subtype_markers,	
				params.major_markers,				
				tier_two.out)
				
			QC_subtype(subtypes_exporter.out,
					params.subtype_method_nostrat,
					'subtypes',
	                params.subtype_method)
		}
	// emit:
	//	subtypes_exporter:.out
}

workflow STRATIFY {

	take: out
		  mask_overlay
		  tier_one
	main:
		cellFiles=get_cell_files(params.deep_imcyto, params.cellprofiler)
		rawMasks=get_imcyto_raw_masks()
		
		if(params.exclude_cell_lineage != 'None') {
			tier_one_ref(
				params.major_method,
				'major',
				params.major_markers,
				params.major_markers,
				true,
				"${params.output_dir}/nfData", 
				out
			)
		} else {
			tier_one_ref(
				params.major_method,
				'major',
				params.major_markers,
				params.major_markers,
				true,
				"${params.output_dir}/nfData", 
				tier_one
			)
		}
		
		//if(params.imcyto) {
		//	csm_submit(
		//		rawMasks.combine(cellFiles, by: 0), 
		//		"${params.output_dir}/major/${params.major_markers}/csm"
		//	)
		//	csm_export(
		//		"${params.output_dir}/major/${params.major_markers}/csm",
		//		csm_submit.out.collect())
		//	build_strata_model(
		//		tier_one, 
		//		tier_one_ref.out,
		//		csm_export.out,
		//		mask_overlay
		//	)
		//} else {
		
			build_strata_model(
				tier_one, 
				tier_one_ref.out,
				out,
				mask_overlay
			)
		//}		 		
	emit:
		build_strata_model.out
}

workflow SUBSAMPLING {

	///////////////////////////////////////////////
	/* --        ROBUSTNESS ANALYSES        -- */
	//////////////////////////////////////////////
	take: out
	main:
		
		if(params.sampled || params.clustered) 	{
			call_cluster(
				params.subtype_method, 
				"major",
				params.annotate_markers,
				params.annotate_markers,
				false,
				"${params.output_dir}/nfData",
				out
			)
		 	major_exporter(
				params.subtype_method,
				'major',
				params.annotate_markers,
				params.annotate_markers,
				call_cluster.out.collect())
				
			QC_major(
				major_exporter.out,
				"${params.annotate_markers}",
				'major',
				params.subtype_method)
		}
		if(params.sampled) {
			pars = iterations.combine(subsample_methods)
					.combine(Chanel.of(["major"]))
					.combine(subsample_markers)
			call_subsampled(
	            pars,
	            "${params.output_dir}/nfData",
				out
			)

			match_clusters(
				params.output_dir, 
				subsample_methods,
				subsample_markers,
				params.annotate_markers,
				call_subsampled.out.collect(), 
				call_cluster.out.collect()
			)
			plot_subsampled(match_clusters.out.collect())
		}

		//plot_dr(
		//	"${params.dimred_method}", 
		//	"${params.output_dir}/nfData", 
		//	call_cluster.out.collect()
		//)
}
