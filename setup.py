from distutils.core import setup
import setuptools



with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
   name='aliasfinder',
   version='1.0',
   description='The AliasFinder is a Python script for uncomplicated '
               'alias testing based on the method introduced by '
               'Dawson & Fabrycky (2010)',
   license="MIT",
   long_description=long_description,
   long_description_content_type="text/markdown",
   author='Stephan Stock, Jonas Kemmer @ ZAH, Landessternwarte Heidelberg',
   author_email='jkemmer@lsw.uni-heidelberg.de',
   url="https://github.com/JonasKemmer/AliasFinder",
   packages=setuptools.find_packages(),
   install_requires=['astropy', 'numpy', 'matplotlib >= 3.1.1',
                     'scipy',
                     'pyyaml', 'tqdm'],
   python_requires='>=3.6',
   entry_points={'console_scripts': ['AliasFinder=aliasfinder:main']},
   )
