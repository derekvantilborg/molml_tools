from setuptools import setup


setup(
    name='molml_tools',
    version='1.0.7',
    packages=['molml', 'molml.Viz', 'molml.Data', 'molml.Data.Chembl', 'molml.Data.Data_prep',
              'molml.Data.example_data', 'molml.Tools', 'molml.Datastructures', 'molml.Representations'],
    url='https://github.com/derekvantilborg/molmlkit',
    license='MIT',
    author='Derek van Tilborg',
    author_email='d.w.v.tilborg@tue.nl',
    description='A collection of tools for molecular machine learning',
    install_requires=[
        'tqdm',
        'requests',
        'twine',
        'importlib-metadata',
        'pandas',
        'numpy',
        'matplotlib',
        'seaborn',
        'chembl_webresource_client',
        'scikit-learn',
        'scikit-optimize',
        'matplotlib',
        'rdkit-pypi',
        'keras'
    ],
    include_package_data=True,
    package_data={'': ['Data/example_data/*']}
)
