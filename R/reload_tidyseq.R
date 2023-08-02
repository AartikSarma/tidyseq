#Function for reloading while debugging 

reload_tidyseq <- function( ){
  detach("package:tidyseq", unload = T)
  library(tidyseq)
}