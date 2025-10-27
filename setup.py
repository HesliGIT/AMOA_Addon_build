from setuptools import setup, Extension
import pybind11
import sys

cpp_args = []
link_args = []

if sys.platform == 'win32':
    cpp_args = [
        '/std:c++17',
        '/O2',
        '/openmp',
    ]
    link_args = ['/openmp']
else:
    cpp_args = [
        '-std=c++17',
        '-O3',
        '-march=native',
        '-fopenmp',
    ]
    link_args = ['-fopenmp']

ext_modules = [
    Extension(
        'amoa_cpp',
        ['spatial_grid.cpp'],
        include_dirs=[pybind11.get_include()],
        language='c++',
        extra_compile_args=cpp_args,
        extra_link_args=link_args,
    ),
]

setup(
    name='amoa_cpp',
    version='1.0.1',
    author='Hesli Reiling',
    description='A C++ accelerated module for AMOA',
    long_description='This package provides a high-performance C++ backend for calculating object similarity, using batching and parallel processing.',
    ext_modules=ext_modules,
    python_requires='>=3.6',
    install_requires=[
        'pybind11>=2.6'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'Programming Language :: C++',
    ],
)