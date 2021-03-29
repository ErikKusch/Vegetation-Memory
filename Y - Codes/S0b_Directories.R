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
Dir.ERA <- paste(Dir.Data, "/2 - ERA5-Land", sep="")
if(!dir.exists(Dir.ERA)){dir.create(Dir.ERA)}
# # WORKING DIRECTORY FOR PROCESSED ERA5 DATA
# Dir.ERA.Monthly <- paste(Dir.Data, "/2b - ERA5_Monthly", sep="")
# if(!dir.exists(Dir.ERA.Monthly)){dir.create(Dir.ERA.Monthly)}
# # WORKING DIRECTORY FOR KRIGING COVARIATES
# Dir.KrigCov <- paste(Dir.Data, "/3 - Kriging", sep="")
# # WORKING DIRECTORY FOR MEMORY EFFECT DATA
Dir.Memory <- paste(Dir.Data, "/4 - Memory_Effects", sep="")
if(!dir.exists(Dir.Memory)){dir.create(Dir.Memory)}
# WORKING DIRECTORY FOR SHAPEFILES (contains masking file for water bodies)
Dir.Mask <- paste(Dir.Data, "/5 - ShapeFiles", sep="")
# WORKING DIRECTORY FOR PLOTS
Dir.Plots <- paste(mainDir, "/Z - Plots", sep="")
if(!dir.exists(Dir.Plots)){dir.create(Dir.Plots)}
