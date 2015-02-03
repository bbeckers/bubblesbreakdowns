# bubblesbreakdowns
This git stores the R-code used to download US realtime data, combine it with unrevised data and make a pseudo-out-of-sample
realtime experiment for the prediction of inflation and industrial production with a time horizon of 12-month.

## getting data
"Getting and preparing realtime data.R" downloads FRED Realtime database and binds monthly data that have a vintage prior to
1998 (exclusive) to sets. Each set contains the remaining nine variables as the have been avaible at each vintage month. 
There are ?? sets saved as sets. Produces "sets.RData", "overview.csv" (containing an overview of the 9 remaining) and "overview.complete.csv" (overview of all monthly variables).
