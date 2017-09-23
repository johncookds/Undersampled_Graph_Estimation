from setuptools import setup

with open('README.md') as description:
    description = description.read()

#execfile('estimation/version.py') # Acquire version constants.

# Define some package entry points. These will be command-line scripts that get
# installed into the user's PATH

install_requires = []

dependency_links=[]

setup(
    name='Undersampled_Graph_Estimation',
    description='Tools to explore dynamic causal graphs in the case of undersampled data',
    version='0.0.1',
    author='John Cook',
    author_email='johncookds@gmail.com',
    long_description=description,
    include_package_data=True, # Include files listed in MANIFEST.in
    packages=['dbnestimation'], # Sub-packages must be explicitly listed.
    #entry_points=epoints,
    install_requires=install_requires, # List of dependencies.
    dependency_links=dependency_links,
    zip_safe=False
    ) # Override annoying default behavior of easy_install.
