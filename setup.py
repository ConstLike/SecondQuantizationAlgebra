#!/usr/bin/env python3
"""
Setup script for Second Quantization Algebra (SQA) package.
"""

from setuptools import setup, find_packages
import os

# Read README for long description
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ''

setup(
    name='sqa',
    version='3.0.0',
    description='Second Quantization Algebra for quantum chemistry',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    author='Eric Neuscamman, Masaaki Saitow, Konstantin Komarov',
    author_email='constlike@gmail.com',
    maintainer='Konstantin Komarov',
    maintainer_email='constlike@gmail.com',
    url='https://github.com/ConstLike/SecondQuantizationAlgebra',
    license='Academic/Research Use',

    # Package discovery
    py_modules=[
        'secondQuantizationAlgebra',
        'sqaOptions',
        'sqaIndex',
        'sqaSymmetry',
        'sqaTensor',
        'sqaTerm',
        'sqaNormalOrder',
        'sqaCommutator',
        'sqaMisc',
        'sqaMisc2',
        'sqaDecomposition_sf',
        'sqaDecomposition_so',
        'sqaConvert',
        'sqaConvert2',
        'sqaSLee',
    ],

    # Python version requirement
    python_requires='>=3.6',

    # Package metadata for PyPI
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: Free for non-commercial use',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],

    keywords=[
        'quantum chemistry',
        'second quantization',
        'symbolic algebra',
        'tensor algebra',
        'coupled cluster',
        'configuration interaction',
        'density matrices',
        'multireference',
        'normal ordering',
    ],

    # Project URLs
    project_urls={
        'Documentation': 'https://github.com/ConstLike/SecondQuantizationAlgebra/blob/master/SQA_User_Guide.md',
        'Source': 'https://github.com/ConstLike/SecondQuantizationAlgebra',
        'Tracker': 'https://github.com/ConstLike/SecondQuantizationAlgebra/issues',
    },
)
