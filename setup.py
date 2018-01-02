from setuptools import setup

FNAME_README = "README.rst"
FNAME_REQUIREMENTS = "requirements.txt"
VERSION = "3.1.2"

readme = open(FNAME_README).read()
requirements = open(FNAME_REQUIREMENTS).read().split()

setup(
    name="cosmopy",
    version=VERSION,
    author="Luke Zoltan Kelley",
    author_email="lkelley@cfa.harvard.edu",
    description=("General, commonly used functions for other projects."),
    license="MIT",
    url="https://github.com/lzkelley/cosmopy",
    download_url="https://github.com/lzkelley/cosmopy/archive/v{}.tar.gz".format(VERSION),
    packages=['cosmopy'],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'cosmo = cosmopy.__main__:main'
        ]
    },
    install_requires=requirements,
    long_description=readme,
    keywords=['utilities', 'physics', 'astronomy', 'cosmology',
              'astrophysics', 'calculator'],
    CLASSIFIERS=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Education",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Utilities"
    ]
)
