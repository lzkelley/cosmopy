cosmopy
=======

.. image:: https://img.shields.io/travis/lzkelley/cosmopy.svg
    :target: https://travis-ci.org/lzkelley/cosmopy?branch=master

.. image:: https://img.shields.io/codecov/c/github/lzkelley/cosmopy/master.svg
    :target: https://codecov.io/gh/lzkelley/cosmopy

Quickly calculate cosmological parameters.

- Provides both command-line and API interfaces.

- Uses the machinery from the `astropy.cosmology` package.

The user provides an input parameter (e.g. redshift or luminosity-distance) and recieves the entire set of cosmological measures corresponding to the epoch thus specified.


The below gif shows three examples: inputting a redshift, a luminosity distance (`-dl 400Mpc`), and an age of the universe (`-ta 3.2Gyr`).

.. image:: https://raw.githubusercontent.com/lzkelley/cosmopy/dev/docs/cosmopy_demo.gif
   :height: 600px


Installation
------------

Using `pip`:

::

    pip install cosmopy

From source:

::

 git clone git@github.com:lzkelley/cosmopy.git
 pip install cosmopy


Usage
-----
`cosmopy` can be used via the command-line `cosmo` command, or via python API by importing the module directly.

- Command Line:

    ::

        $ cosmo --help

         usage: cosmo [-h] [-z Z] [-a A] [-dc DC] [-dl DL] [-tl TL] [-ta TA] [-v]

         cosmopy: cosmological calculator.

         optional arguments:
           -h, --help      show this help message and exit
           -z Z            target redshift z
           -a A            target scale factor a
           -dc DC, -cd DC  target coming distance D_C
           -dl DL, -ld DL  target luminosity distance D_L
           -tl TL, -lt TL  target look-back time T_L
           -ta TA, -at TA  target universe age T_A
           -v, --version   print version information.

    For example, entering an input redshift of 0.2:

    ::

        $ cosmo -z 0.2

               z = 0.2000                                : Redshift
               a = 0.8333                                : Scale-factor
             D_c = 815.3960 Mpc      ~ 2.5160e+27 cm     : Comoving Distance
             D_L = 978.4752 Mpc      ~ 3.0193e+27 cm     : Luminosity Distance
             D_A = 679.4967 Mpc      ~ 2.0967e+27 cm     : Angular Diameter Distance
          Arcsec = 3294.2928 pc      ~ 1.0165e+22 cm     : Arcsecond Scale
            T_lb = 2.4277 Gyr        ~ 7.6613e+16 s      : Lookback Time
             T_a = 11.3235 Gyr       ~ 3.5734e+17 s      : Age of the Universe
              DM = 39.9527                               : Distance Modulus

    or an input luminosity-distance of 400 Mpc:

    ::

        $ cosmo -dl 400Mpc

               z = 0.0880                                : Redshift
               a = 0.9192                                : Scale-factor
             D_c = 367.6631 Mpc      ~ 1.1345e+27 cm     : Comoving Distance
             D_L = 400.0000 Mpc      ~ 1.2343e+27 cm     : Luminosity Distance
             D_A = 337.9403 Mpc      ~ 1.0428e+27 cm     : Angular Diameter Distance
          Arcsec = 1638.3809 pc      ~ 5.0555e+21 cm     : Arcsecond Scale
            T_lb = 1.1496 Gyr        ~ 3.6280e+16 s      : Lookback Time
             T_a = 12.6016 Gyr       ~ 3.9768e+17 s      : Age of the Universe
              DM = 38.0103                               : Distance Modulus

- Python API

    The module can be imported as `cosmopy`, from which the primary access point is the `api` function which accepts two arguments: a `key` (a target cosmological parameter) and a `value` of that parameter (optionally including units).  The function returns a dictionary with the computed values as key: value pairs (both strings).  For example:

    ::

        $ python -c "import cosmopy; print(cosmopy.api('dl', '1.2 Gpc'))"
         {'z': '0.2396', 'dl': '1200.0000 Mpc', 'tl': '2.8359 Gyr', 'dc': '968.0336 Mpc', 'ta': '10.9153 Gyr', 'da': '780.9075 Mpc', 'dm': '40.3959', 'arc': '3785.9464 pc', 'a': '0.8067'}
