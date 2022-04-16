# setup.py: The setup file
from setuptools import setup, Extension


setup(name='mykmeanssp',
      version='1.1',
      author="Nader, Weam",
      description='kmeans algorithem for sp class',
      ext_modules=[Extension('mykmeanssp',['main.c',  'spkmeansmodule.c'])])
