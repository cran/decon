.onLoad <- function(lib, pkg)
   packageStartupMessage("Decon 1.2-1 loaded")

.onUnload <- function(libpath)
    library.dynam.unload("decon",  libpath)

