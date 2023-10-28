include { qc_select_images; qc_create_single_channel_images; 
		qc_overlay; qc_intensity_heatmap} from '../modules/qc.nf'
		
workflow QC {

	take: out
		  markers
		  subset
		  method
	main:
	
			run_dirs = Channel.fromPath(params.image_dir, type: "dir", relative: false)
			run_dirs.view()
				qc_select_images(
					markers, 
					subset,
					method,
					out
				)
				qc_create_single_channel_images(
					markers,
					subset,
					method,
					run_dirs,
					qc_select_images.out.collect()
				)
				run_dirs = Channel.fromPath(params.input_dir, type: "dir", relative: false)
				qc_overlay(
					markers,
					subset, 
					method,
					run_dirs,
					qc_create_single_channel_images.out.collect()
			)
		qc_intensity_heatmap(
			markers,
			subset, 
			method,
			out
		)
}
