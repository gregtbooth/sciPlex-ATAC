library(stringr)
library(Matrix)

options(stringsAsFactors = FALSE)

write_vector = function(vec, file) {
	fileConn<-file(file)
	writeLines(vec, fileConn)
	close(fileConn)
}

command_exists = function(command) {
	return( 0 == system(paste0('which ', command), ignore.stdout=TRUE, ignore.stderr=TRUE))
}

write_mtx_file = function(mat, filename) {
	if (str_detect(filename, '[.]gz$')) {
		uncompressed_file = str_replace(filename, '[.]gz$', '')
		Matrix::writeMM(mat, file=uncompressed_file)

		if (command_exists('pigz')) {
			system(paste0('pigz -f ', uncompressed_file))
		} else {
			system(paste0('gzip ', uncompressed_file))
		}
	} else {
		# Just write the uncompressed file
		Matrix::writeMM(mat, file=filename)
	}

	# Write files for rows/cols
	dim_files = .get_aux_files(filename)

	write_vector(rownames(mat), dim_files$features)
	write_vector(colnames(mat), dim_files$cells)
}

.get_aux_files = function(mtx_file) {
	base_name = str_replace(mtx_file, '[.]gz$', '')
	base_name = str_replace(base_name, '[.]mtx$', '')

	features_file = paste0(base_name, '.rows.txt')
	cells_file = paste0(base_name, '.columns.txt')

	return(list(features=features_file, cells=cells_file))
}

load_mtx_file = function(mtx_file) {
	mat = readMM(mtx_file)
	dim_files = .get_aux_files(mtx_file)

	if (! file.exists(dim_files$features)) {
		stop(paste0(dim_files$features, ' file not found when loading ', mtx_file))
	}

	if (! file.exists(dim_files$cells)) {
		stop(paste0(dim_files$cells, ' file not found when loading ', mtx_file))
	}

	rownames(mat) = read.delim(dim_files$features, header=FALSE)$V1
	colnames(mat) = read.delim(dim_files$cells, header=FALSE)$V1
	return(mat)
}