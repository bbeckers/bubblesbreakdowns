This folder contains three programs:

Two one-sided HP filters:
one_sided_hp_filter_kalman.m
one_sided_hp_filter_serial.m

And one standard (two-sided) HP-Filter
hp_filter_sparse.m

The standard (two-sided) HP-Filter was a byproduct of one_sided_hp_filter_serial.m and 
provides a fast implementation of the standard HP filter using sparse matrices that might 
be useful (especially iterative pricedures).

The two one-sided HP filters yield numerically identical results (at least during my 
testing). While the Kalman filter implementation (one_sided_hp_filter_kalman.m) is
perhaps more appealing from a filtering perspective, the serial implementation of the 
two-sided filter (one_sided_hp_filter_serial.m) is sveral times faster.

This software is provided free of charge with no warranty or guarranty that it will work,
but with the hope that it will provide useful. You are free to use/modify/redistribute these 
programe so long as original authorship credit is given and you in no way impinge on its 
free distribution.

Alexander Meyer-Gohde
Humboldt-University at Berlin
Sept. 2010