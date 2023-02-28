install.packages('neonUtilities')

library(neonUtilities)
library(tidyverse)

secchi <- loadByProduct(dpID = 'DP1.20252.001',
              site = c('BARC', 'CRAM', 'LIRO', 'PRLA', 'PRPO', 'SUGG', 'TOOK'), 
              startdate = '2021-01', enddate = '2022-12', package = 'basic')

dep_secchi <- secchi$dep_secchi |> 
  select(siteID, date, clearToBottom, secchiMeanDepth, maxDepth) |> 
  mutate(secchiMeanDepth = ifelse(clearToBottom == 'Y',
                                  maxDepth, secchiMeanDepth))

ggplot(dep_secchi, aes(x=date, y=secchiMeanDepth)) +
  geom_point()+
  facet_wrap(~siteID, scales = 'free')

dep_secchi |> 
  group_by(siteID) |> 
  summarise(max_Sd = max(secchiMeanDepth, na.rm = T),
            min_Sd = min(secchiMeanDepth, na.rm = T)) |> 
  mutate(kmax = 1.7/max_Sd,
         kmin = 1.7/min_Sd)


dep_secchi |> 
  group_by(siteID) |> 
  summarise(mean_Sd = mean(secchiMeanDepth, na.rm = T)) |> 
  mutate(k = 1.7 / mean_Sd)

