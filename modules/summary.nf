
process compare_dp {

	 tag "sub"
        label 'medium_mem'
        maxRetries 1

        input:
                val method
                val exporter

        script:
        """
                export BASE_DIR=$baseDir
        """

}

process compare_subsampled {

	tag "sub"
        label 'xs'
        maxRetries 1

        input:
                val method
                val exporter

        script:
        """
                export BASE_DIR=$baseDir
        """
}

process compare_tcra {
		tag "sub"
        label 'xs'
        maxRetries 1

        input:
                val method
                val exporter

        script:
        """
                export BASE_DIR=$baseDir
        """
}

process compare_flow {
        tag "sub"
        label 'xs'
        maxRetries 1

        input:
                val method
                val exporter

        script:
        """
                export BASE_DIR=$baseDir
        """
}

process compare_panels {

	tag "sub"
        label 'xs'
        maxRetries 1

        input:
                val method
                val exporter

        script:
        """
                export BASE_DIR=$baseDir
        """

}
