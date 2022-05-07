from setuptools import setup, find_packages

setup(
    name='netprop',
    version='0.0.5.1',
    packages=find_packages(),
    install_requires=[
        'matplotlib==3.5.1',
        'mygene==3.2.2',
        'ndex2==3.4.0',
        'networkx==2.6.3',
        'numpy==1.21.4',
        'pandas==1.3.4',
        'pydantic==1.8.2',
        'scipy==1.7.3',
        'sortedcontainers==2.4.0'
    ]
)