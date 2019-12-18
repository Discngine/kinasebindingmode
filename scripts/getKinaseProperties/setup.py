from setuptools import setup, find_packages

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
        name='getKinaseProperties',
        version='0.1.0',
        description='Python script for kinase structural features processing, in order to classify potential binding modes',
        lon_description=readme,
        author='Anne-Sophie',
        author_email='',
        url='https://github.com/Discngine/kinasebindingmode/getKinaseProperties',
        license=license,
        package=find_packages(exclude=('tests', 'docs'))
)
