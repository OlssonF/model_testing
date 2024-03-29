lake <- 'fcre'
model <- 'Simstrat'

setwd(file.path(here::here(), lake))

LakeEnsemblR::export_meteo(config_file = paste0('LakeEnsemblR_', lake,'.yaml'), model = model)
LakeEnsemblR::export_extinction(config_file = paste0('LakeEnsemblR_', lake,'.yaml'), model = model)

LakeEnsemblR::export_config(config_file = paste0('LakeEnsemblR_', lake,'.yaml'), model = model)
