from setuptools import setup

readme = open('README.rst').read()

requirements = []

setup(
    name="cosmocalc",
    version=2.0,
    author="Luke Zoltan Kelley",
    author_email="lkelley@cfa.harvard.edu",
    description=("General, commonly used functions for other projects."),
    license="MIT",
    keywords="",
    url="https://bitbucket.org/lzkelley/zcode/",
    packages=['cosmocalc'],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'cosmo = cosmocalc.__main__:main'
        ]
    },
    install_requires=requirements,
    long_description=readme,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
    ],
)
