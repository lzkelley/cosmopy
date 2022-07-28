from setuptools import setup

FNAME_README = "README.rst"
FNAME_REQUIREMENTS = "requirements.txt"
FNAME_VERSION = "cosmopy/VERSION.txt"

with open(FNAME_REQUIREMENTS, "r") as inn:
    requirements = inn.read().splitlines()

with open(FNAME_README, "r") as inn:
    long_description = inn.read().strip()

with open(FNAME_VERSION, "r") as inn:
    version = inn.read().strip()

setup(
    name="cosmopy",
    author="Luke Zoltan Kelley",
    author_email="lzkelley@northwestern.edu",
    version=version,
    description="",
    license="MIT",
    url="https://github.com/lzkelley/cosmopy",
    download_url="https://github.com/lzkelley/cosmopy/archive/v{}.tar.gz".format(version),
    packages=['cosmopy'],
    entry_points={
        'console_scripts': [
            'cosmo = cosmopy.__main__:main'
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords=['utilities', 'physics', 'astronomy', 'cosmology', 'astrophysics', 'calculator'],
    python_requires=">=3.5",
)
