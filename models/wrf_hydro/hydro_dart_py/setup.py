from setuptools import find_packages, setup

setup(
    name='hydrodartpy',
    version='0.0.1',
    packages=find_packages(),
    package_data={'hydrodartpy': ['core/data/*']},
    url='https://github.com/NCAR/wrf_hydro_dart',
    license='MIT',
    install_requires=[
        'boltons',
        'datetime',
        'deepdiff',
        'f90nml',
        'netCDF4',
        'pandas',
        'pathlib',
        'pytest',
        'pytest-html',
        'pytest-datadir-ng',
        'pywrfhydro',
        'wrfhydropy',
        'xarray',
        'ruamel.yaml'
    ],
    author='James McCreight',
    author_email='jamesmcc@ucar.edu',
    description='API for wrf_hydro_dart',
)
