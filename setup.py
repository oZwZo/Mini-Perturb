import os
from setuptools import setup
from setuptools import find_packages

setup(
    name = 'MiniPert',
    version = '0.1.1',
    description = 'MiniPert: Infering perturbation effect for multi-batch mini-bulk RNA-seq data',
    author = 'Weizhong Zheng',
    author_email = "zhengwzh@connect.hku.hk",
    packages = ['MiniPert'],
    install_requires = [
        'numpy',
        'pandas',
        'scanpy>=1.6.0',
        'matplotlib',
        'seaborn',
        'scipy',
        'statsmodels',
        'dtw',
        'scikit-fda==0.7.1',
        'tqdm'
    ],
    
)