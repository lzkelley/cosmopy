cosmopy
=========

.. image:: https://travis-ci.org/lzkelley/cosmopy.svg?branch=master
    :target: https://travis-ci.org/lzkelley/cosmopy?branch=master
    
.. image:: https://coveralls.io/repos/github/lzkelley/cosmopy/badge.svg?branch=master
    :target: https://coveralls.io/github/lzkelley/cosmopy?branch=master
    
Quickly calculate cosmological parameters.  

- Provides both command-line and API interfaces.  

- Uses the machinery from the `astropy.cosmology` package.

The user provides an input parameter (e.g. redshift or luminosity-distance) and recieves the entire set of cosmological measures corresponding to the epoch thus specified.

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

::
    
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
