include { qc_select_images; qc_create_single_channel_images; 
		qc_overlay; qc_intensity_heatmap} from '../modules/qc.nf'
		
workflow QC {

	take: out
		  markers
		  subset
		  method
	main:
	
		//if(file(params.image_dir).isDirectory())	{
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
					qc_select_images.out
				)
				qc_overlay(
					markers,
					subset, 
					method,
					qc_create_single_channel_images.out
			)
		//}
		qc_intensity_heatmap(
			markers,
			subset, 
			method,
			out
		)
}
