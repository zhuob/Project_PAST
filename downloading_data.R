### Downloading Arabidopsis data from NCBI
getwd()

download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE56922&format=file&file=GSE56922%5Fraw%5Fcounting%5Fdata10%2Ecsv%2Egz", 
              destfile = "GSE56922_raw_counting_data10.csv.gz")


untar("GSE56922_raw_counting_data10.csv.gz", list = TRUE)
