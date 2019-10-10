mainDir <- getwd() # extract the project folder location
# WORKING DIRECTORY FOR CODES
Dir.Codes <- paste(mainDir, "/Y - Codes", sep="")
# WORKING DIRECTORY FOR DATA
Dir.Data <- paste(mainDir, "/X - Data", sep="")
# WORKING DIRECTORY FOR RAW GIMMS DATA
Dir.Gimms <- paste(Dir.Data, "/1 - GIMMs_Raw", sep="")
if(!dir.exists(Dir.Gimms)){dir.create(Dir.Gimms)}
# WORKING DIRECTORY FOR PROCESSED GIMMS DATA
Dir.Gimms.Monthly <- paste(Dir.Data, "/1b - GIMMs_Monthly", sep="")
if(!dir.exists(Dir.Gimms.Monthly)){dir.create(Dir.Gimms.Monthly)}
# WORKING DIRECTORY FOR RAW ERA5 DATA
Dir.ERA <- paste(Dir.Data, "/2 - ERA5_Raw", sep="")
# WORKING DIRECTORY FOR PROCESSED ERA5 DATA
Dir.ERA.Monthly <- paste(Dir.Data, "/2b - ERA5_Monthly", sep="")
if(!dir.exists(Dir.ERA.Monthly)){dir.create(Dir.ERA.Monthly)}
# WORKING DIRECTORY FOR KRIGING COVARIATES
Dir.KrigCov <- paste(Dir.Data, "/3 - Kriging", sep="")
# WORKING DIRECTORY FOR MEMORY EFFECT DATA
Dir.Memory <- paste(Dir.Data, "/4 - Memory_Effects", sep="")
# WORKING DIRECTORY FOR COMPADRE DATA
Dir.Compadre <- paste(Dir.Data, "/5 - COMPADRE", sep="")
# WORKING DIRECTORY FOR TRY PFT DATA
if(!dir.exists(Dir.Memory)){dir.create(Dir.Memory)}
Dir.TRY <- paste(Dir.Data, "/5 - TRY", sep="")
# WORKING DIRECTORY FOR OCCURENCE DATA
Dir.OCCs <- paste(Dir.Data, "/6 - Occurences", sep="")
if(!dir.exists(Dir.OCCs)){dir.create(Dir.OCCs)}
# WORKING DIRECTORY FOR SHAPEFILES (contains masking file for water bodies)
Dir.Mask <- paste(Dir.Data, "/7 - ShapeFiles", sep="")