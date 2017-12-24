"""
"""


class Settings:
    """Object to store the current configuration of CosmoCalc.
    """

    # Print excess output to stdout
    verbose = True
    # Print the cosmological parameters
    print_flag = False
    # Default whether to rebuild integration table
    build_flag = False

    # Use input distance in parsecs
    use_pc = False
    # Use input distance in Megaparsecs
    use_mpc = False
    # Use input distance in lightyears
    use_ly = False

    # Use input time     in years
    use_yr = False
    # Use input time     in megayears
    use_myr = False

    # Target epochs (calculate epoch based on given parameter below)
    # Redshift (z)
    z = None
    # Scale Factor (a)
    a = None
    # Comoving Distance
    cd = None
    # Luminosity Distance
    ld = None
    # Lookback time
    lt = None
    # Age of the universe (time)
    t = None

    # Root finding tolerance
    root_tol = 1.0e-10
    # Number of iterations for quadrature
    quad_iter = 1000
