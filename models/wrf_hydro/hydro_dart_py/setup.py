from setuptools import find_packages, setup

setup(
    name='hydrodartpy',
    version='0.0.1',
    packages=find_packages(),
    package_data={'hydrodartpy': ['core/data/*']},
    url='https://github.com/NCAR/DART',
    license='MIT',
    install_requires=[
        'wrfhydropy>=0.0.18',
        'ruamel.yaml'
    ],
    author='James McCreight',
    author_email='jamesmcc@ucar.edu',
    description='API for wrf_hydro_dart',
)
