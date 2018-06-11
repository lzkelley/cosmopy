CHANGES
=======


To-Do
-----
    -   Improve tests to better capture install process (requirements, etc).
    -   CMB Temperature
    -   Allow lists of target values (instead of single).
    -   Online documentation.
    

Current
-------
    -   Added redshift to comoving and luminosity distance methods.
    -   Added differential volume calculation (dVc/dz)


v3.2 - 2018-01-03
-------------------


v3.1.3 - 2018-01-03
-------------------
    -   Added colored output using the 'click' package (added to dependencies).
    -   Added significant unit-testing.  Coverage currently at 98%.
    -   Added function documentation.
    

v3.1.2 - 2018-01-01
-------------------
    -   Fix issue with version file.  Use hardcoded version strings in both `setup.py` and `__init__.py` for now.
    -   Fix python2-python3 incompatibility issue regarding use of `super`.
        -   Added `future` package to requirements to deal with this.
    -   Improve the README usage documentation.


v3.1.1 - 2018-01-01
-------------------
    -   Fix install issues from missing 'requirements' file (wasn't being loaded properly in 'setup.py').
    -   When no arguments are given, explain that arguments are required and print help message.
    -   Added version (and other metadata) info to '__init__.py'.


v3.1 - 2017-12-31
-----------------
    -   Renamed everything from `cosmocalc` to `cosmopy` (former already taken on pypi).
    -   Added travis and coveralls integration and badge-icons.
    -   Command-line input is parsed using astropy to identify input units.
    -   API functions added to retrieve cosmological parameters via python method calls.
        -   Can also retrieve underlying `Cosmology` instances directly.
    

v3.0 - 2017-12-24
-----------------
    - Completely restructured and rewrote code to use an `astropy.cosmology.FlatLambdaCDM` subclass for bulk of cosmological calculations.  This class uses stored grids of distance measures to perform interpolation when inverting from a distance measure to a redshift.
    - The basic functionality is working in which an input redshift, scale-factor, distance (lum or com), of time (age or lookback) can be input, and all of the other parameters will be calculated and printed.


    - cosmopy/
        - __main__.py
            - This provides the entire API for command-line accessible functions (at the moment).
            - `parse_input()` [NEW-FUNCTION]
                - Take an input string which specifies a quantitative value and convert it into a usable numerical object using `astropy.units.Quantity`.
            - `calc()`
                - Runs all of the basic cosmological calculations.
            - `output()`
                - Format and print the resulting values.
        - cosmology.py [NEW-FILE]
            - Cosmology [NEW-CLASS]
                - Currently uses fixed set of cosmological parameters (i.e. Omega values) and calculated distance measures.
        - tests/ [NEW-FOLDER]
            - test_main.py [NEW-FILE]
                - Added a few very simple tests for parsing input values.

    - parameters.py [DELETED]
    - settings.py [DELETED]


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
