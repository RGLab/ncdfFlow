# LdFlags <- function(){
#   os.type <- Sys.info()[['sysname']]
#   libpath <- paste0("lib", Sys.getenv("R_ARCH"), "/libncdfFlow.a")
# 
#   flags <- tools::file_path_as_absolute(base::system.file(libpath, package = "ncdfFlow" ))
#   #add flags for hdf5 lib
#   if(os.type != "Windows")
#   {
#     dll <- system.file("libs/ncdfFlow.so", package = "ncdfFlow")
#     hdf.paths <- system2("ldd" , dll, stdout = TRUE)
#     hdf.paths <- hdf.paths[grep("libhdf5", hdf.paths)]
#     if(os.type == "Linux")
#     {
#       hdf.paths <- strsplit(hdf.paths, " ")[[1]][3]  
#     }else
#       hdf.paths <- strsplit(hdf.paths, " ")[[1]][1]  
#     
#     hdf.paths <- dirname(hdf.paths)
#     flags <- paste0(flags, " -L", hdf.paths, " -lhdf5")
#   }
#     
#   cat(flags)
#   
# }
