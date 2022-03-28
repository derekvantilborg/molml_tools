from setuptools import setup

setup(
    name='Molml-tools',
    version='1.0.0',
    packages=['Viz', 'Data', 'Data.Chembl', 'Data.Data_prep', 'Tools', 'example_data', 'Representations'],
    url='https://github.com/derekvantilborg/Molml-tools',
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
        'rdkit-pypi'
    ],
    include_package_data=True,
    package_data={'': ['example_data/*']}
)
