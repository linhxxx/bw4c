from setuptools import setup
from Cython.Build import cythonize

setup(
    name='bw4c Modules',
    ext_modules=cythonize("bw4cModule.pyx"),
    zip_safe=False,
)
