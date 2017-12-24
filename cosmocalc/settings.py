"""
"""


class Settings:
    """Object to store the current configuration of CosmoCalc.
    """

    verbose = True   # Print excess output to stdout
    print_flag = False  # Print the cosmological parameters
    build_flag = False  # Default whether to rebuild integration table

    use_pc = False  # Use input distance in parsecs
    use_mpc = False  # Use input distance in Megaparsecs
    use_ly = False  # Use input distance in lightyears

    use_yr = False  # Use input time     in years
    use_myr = False  # Use input time     in megayears

    # Target epochs (calculate epoch based on given parameter below)
    z = None          # Redshift (z)
    a = None          # Scale Factor (a)
    cd = None          # Comoving Distance
    ld = None          # Luminosity Distance
    lt = None          # Lookback time
    t = None          # Age of the universe (time)

    root_tol = 1.0e-10   # Root finding tolerance
    quad_iter = 1000      # Number of iterations for quadrature
