"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path, listdir
from sys import exit, version

here = path.join(path.abspath(path.dirname(__file__)), 'vairo')

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    # This is the name of your project. The first time you publish this
    # package, this name will be registered for you. It will determine how
    # users can install this project, e.g.:
    #
    # $ pip install sampleproject
    #
    # And where it will live on PyPI: https://pypi.org/project/sampleproject/
    #
    # There are some restrictions on what makes a valid project name
    # specification here:
    # https://packaging.python.org/specifications/core-metadata/#name
    name='vairo',  # Required

    # Versions should comply with PEP 440:
    # https://www.python.org/dev/peps/pep-0440/
    #
    # For a discussion on single-sourcing the version across setup.py and the
    # project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='1.0.0',  # Required

    # This is a one-line description or tagline of what your project does. This
    # corresponds to the "Summary" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#summary
    description='VAIRO guides predictions towards particular dynamic states selecting the prior information input or analysing the results of the search.',
    # Required

    # This is an optional longer description of your project that represents
    # the body of text which users will see when they visit PyPI.
    #
    # Often, this is the same as your README, so you can just read it in from
    # that file directly (as we have already done above)
    #
    # This field corresponds to the "Description" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#description-optional
    long_description=long_description,  # Optional

    # This should be a valid link to your project's main homepage.

    # This field corresponds to the "Home-Page" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#home-page-optional
    url='http://chango.ibmb.csic.es',  # Optional

    # This should be your name or the name of the organization which owns the
    # project.
    author='Isabel Uson',  # Optional

    # This should be a valid email address corresponding to the author listed
    # above.
    author_email='bugs-borges@ibmb.csic.es',  # Optional

    # Classifiers help users find your project by categorizing it.
    #
    # For a list of valid classifiers, see
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[  # Optional
        # Operative system info
        'Operating System :: Unix',

        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research ',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish
        'License :: OSI Approved :: BSD License',

        # Specify the Python versions you support here. In particular, ensure
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],

    # This field adds keywords for your project which will appear on the
    # project page. What does your project relate to?
    #
    # Note that this is a string of words separated by whitespace, not a list.
    keywords='crystallography macromolecular',  # Optional

    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    #   py_modules=["my_module"],
    #
    packages=['vairo', 'vairo.ALEPH', 'vairo.ALEPH.aleph',
              'vairo.ALEPH.aleph.core', 'vairo.app'],
    include_package_data=True,
    package_data={  # Optional
        'vairo': ['binaries/*', 'libs/*', 'templates/*', 'README.md'],
    },
    entry_points={  # Optional
        'console_scripts': [
            'VAIRO=vairo.run_vairo:main',
            'vairo=vairo.run_vairo:main',
            'VAIRO-GUI=vairo.app:main',
            'vairo-gui=vairo.app:main'
        ],
    },
    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[
        'absl-py==1.0.0',
        'biopython==1.79',
        'chex==0.0.7',
        'dm-haiku==0.0.9',
        'dm-tree==0.1.6',
        'immutabledict==2.0.0',
        'ml-collections==0.1.0',
        'numpy==1.21.6',
        'scipy==1.7.0',
        'protobuf==3.20.1',
        'pandas==1.3.4',
        'tensorflow==2.9.0',
        'tensorflow-cpu==2.9.0',
        'matplotlib==3.6.2',
        'python-igraph==0.9.10',
        'pyyaml',
        'future',
        'csb',
        'psutil',
        'paramiko',
        'scikit-learn',
        'pickle5',
        'jinja2',
        'flask',
    ],
    python_requires='>=3.6')
