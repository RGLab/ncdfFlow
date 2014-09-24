context("IO test")

test_that("save_ncfs/load_ncfs", {
      nc <- ncdfFlowSet(GvHD[1:2])
      
      output <- tempfile(pattern = "ncfs")
      
      #save to a new folder
      expect_message(save_ncfs(nc, path = output), "Done")
      expect_message(nc1 <- load_ncfs(output), "Done")
      is_equal_flowSet(nc, nc1)
      
      #overwrite exsiting folder
      expect_error(save_ncfs(nc, path = output), "already exists")
      expect_message(save_ncfs(nc, path = output, overwrite = T), "Done")
      
      #save to the non-existing folder
      suppressWarnings(expect_error(save_ncfs(nc, path = "/faked/folder")))
      
      rdsFile <- list.files(output, pattern = ".rds")
      cdfFile <- list.files(output, pattern = ".nc")
      
      #invalid folder
      newRDS <- "tmp.rds"
      file.rename(file.path(output, rdsFile), file.path(output, newRDS))
      expect_error(save_ncfs(nc, path = output, overwrite = T), "doesn't match")
      
      newRDS1 <- "tmp.dd"
      file.rename(file.path(output, newRDS), file.path(output, newRDS1))
      expect_error(save_ncfs(nc, path = output, overwrite = T), "Not a valid")
      
      #restore rds
      file.rename(file.path(output, newRDS1), file.path(output, rdsFile))
      
      newCDF <- "tmp.nc"
      file.copy(file.path(output, cdfFile), file.path(output, newCDF))
      expect_error(save_ncfs(nc, path = output, overwrite = T), "Not a valid")
      
      #remove the redundant nc file
      file.remove(file.path(output, newCDF))
      
      #remove rds file
      file.remove(file.path(output, rdsFile))
      expect_error(load_ncfs(output), "rds file missing")
      #restore rds file
      saveRDS(nc, file.path(output, rdsFile))
      #remove cdf file
      file.remove(file.path(output, cdfFile))
      #restore it
      file.copy(nc@file, file.path(output, cdfFile))
      #nc.trans
      newCDF <- "tmp.nc.trans"
      file.rename(file.path(output, cdfFile), file.path(output, newCDF))
      expect_message(load_ncfs(output), "Done")
            
})