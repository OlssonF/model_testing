lake <- 'PRPO'

setwd(file.path(here::here(), lake))

LakeEnsemblR::export_meteo(config_file = paste0('LakeEnsemblR_', lake,'.yaml'), model = 'Simstrat')
LakeEnsemblR::export_extinction(config_file = paste0('LakeEnsemblR_', lake,'.yaml'), model = 'Simstrat')
