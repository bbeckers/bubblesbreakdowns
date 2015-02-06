# bubblesbreakdowns
This git stores the R-code used to download US realtime data, combines it with unrevised data and make a pseudo-out-of-sample
realtime experiment for the prediction of inflation and industrial production with a time horizon of 12-month.

## getting data
"Getting and preparing realtime data.R" downloads FRED Realtime database and binds monthly data that have a vintage prior to
1998 (exclusive) to sets. Each set contains the remaining nine variables as the have been avaible at each vintage month. 
There are ?? sets saved as sets. Produces "sets.RData", "overview.csv" (containing an overview of the 9 remaining) and "overview.complete.csv" (overview of all monthly variables).

## starting realtime experiment
"realtime experiment.Rmd" loads unrevised data in "unrevised data.csv" and its description in "NonrevData overview.csv" as well as realtime sets in "sets.RData" and the overview of the realtime data "overview"
- it ensures that the row.names, that is, the dates of the observations correspond
- starts an loop over the vintages (377), merging realtime sets with corresponding unrevised data
- the unrevised data are lagged according to their publication lags given in NonrevData overview
- "bmafo.R" is used to identify the best lag structure of an autoregressive model including the target (core inflation rate) and one additional regressors. bmafo.R is described in "olsbmalag.Rmd". Then a 13 month horizon forecast is made and stored in "forecast.all"

## computing confidence bands of the surprise losses
forecast losse.R