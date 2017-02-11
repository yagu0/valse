for (file in list.files("generate_test_data",full.names=TRUE,recursive=TRUE))
{
	file_name_length = nchar(file)
	if (substr(file, file_name_length, file_name_length) == "R")
		source(file)
}
