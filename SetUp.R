#Set up the path for general scripts
GenScriptPath = "Rscripts/"

#Set up the path for project specific scripts
if(!"ProjectScripts" %in% list.dirs(full.names = FALSE)){
  dir.create("ProjectScripts")
}
ProjScriptPath = paste0(getwd(),"/ProjectScripts/")

#Set up the path for results directory
if(!"GeneralResults" %in% list.dirs(full.names = FALSE)){
  dir.create("GeneralResults")
}
GeneralResultsPath = paste0(getwd(), "/GeneralResults/")

# Source required functions
source(paste0(GenScriptPath,"generalFunc.R"))
source(paste0(ProjScriptPath,"projectFunc.R"))

packageF("sva")
packageF("gplots")
packageF("scales")
packageF("devtools")
packageF("magrittr")
select = dplyr::select
filter = dplyr::filter
mutate = dplyr::mutate
