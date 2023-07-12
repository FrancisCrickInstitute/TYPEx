///////////////////////////////////////////////
/* --        INCLUDE SUBWORKFLOWS        -- */
//////////////////////////////////////////////

		
/////////////////////////////////////////////
/* -- LOAD MODULES and FUNCTIONS -- */
include { format_input; collate_features } from '../modules/format.nf'
include { csm_submit; csm_export } from '../modules/csm.nf'
include { TYPE as tier_one; TYPE as tier_two; TYPE as tier_one_ref; 
		  TYPE as tier_two_nostrat; build_strata_model }  from '../modules/typing.nf'
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
				//['Rphenograph']
				//['flowSOM']
			  )
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
		cellFiles=get_cell_files(params.deep_imcyto, params.mccs)
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
					subtype_methods,
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
                    subtype_methods,
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
			subtypes_exporter(subtype_methods, 
				'subtypes',
				params.subtype_markers,
				params.major_markers,
				tier_two.out)
		} else {
			subtype_methods_nostrat = subtype_methods.map{ method -> [ "${method[0]}_FALSE" ]}
			subtypes_exporter(
				subtype_methods_nostrat,
				'subtypes',
				params.subtype_markers,	
				params.major_markers,				
				tier_two.out)
		}
	// emit:
	//	subtypes_exporter:.out
}

workflow STRATIFY {

	take: out
		  mask_overlay
		  tier_one
	main:
		cellFiles=get_cell_files(params.deep_imcyto, params.mccs)
		rawMasks=get_imcyto_raw_masks()
		
		if(params.most_freq_celltype != 'None') {
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
		pars=iterations.combine(subsample_methods)
				.combine(subsample_markers)
		pars.view()
		if(params.sampled || params.clustered) 	{
			call_cluster(
				subtype_methods, 
				"major",
				params.annotate_markers,
				params.annotate_markers,
				false,
				params.stratify_label,
				"${params.output_dir}/nfData",
				out
			)
		 	major_exporter(
				subtype_methods,
				'major',
				params.annotate_markers,
				params.annotate_markers,
				call_cluster.out.collect())
		}
		if(params.sampled) {
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

workflow SUBSAMPLING {
	take: out
	qc_select_images(
		params.major_markers
	)
	qc_create_single_channel_images(
			qc_select_images.out
	)
}
