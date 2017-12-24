CHANGES
=======

To-Do
-----
    -   'Solve()' and or 'CosmologicalParameters()' should really return a dictionary
    -   CMB Temperature
    -   Produce plots of distances over redshift (and time)
    -   Allow negative target values (i.e. look to future)
    -   Add initial guess for redshift (etc)
    -   Determine values between pairs of redshifts (instead of always relative to zero)
    -   Allow lists of target values (instead of single) --- this can be done via API (i.e. import) but not from the command line yet.


Current
-------



v2.0 - 2014-11-11
-----------------
    -   Module can be imported in (i)python and arrays of parameters can be input and solved for.

v1.3 - 2014-01-19
-----------------
    -   Cosmological parameters can be changed via command-line arguments.

v1.2 - 2013-12-27
-----------------
    -   Determines epoch based on any standard cosmological parameter, and returns the other parameters.

v1.1 - 2013-11-03
-----------------
    -   Removed dependence on precalculated tables --- dynamically integrates to calculate parameters for a given redshift (only).

v1.0 - 2013-10-05
-----------------
    -   Creates a table of times and distances over cosmological history, uses table to interpolate to target redshift (only) - and prints parameters then.
