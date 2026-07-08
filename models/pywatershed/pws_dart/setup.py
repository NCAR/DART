from setuptools import find_packages, setup

setup(
    name='pwsdart',
    version='0.0.1',
    packages=find_packages(),
    package_data={'pwsdart': ['core/data/*']},
    url='https://github.com/NCAR/DART',
    license='MIT',
    install_requires=[
        'ruamel.yaml'
    ],
    author='Ishita Srivastava',
    author_email='ishitas@ucar.edu',
    description='API for py_water_dart',
)
