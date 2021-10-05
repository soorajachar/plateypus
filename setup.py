#! /usr/bin/env python3
from setuptools import setup,find_packages
    
install_requires = [
'pandas',
'seaborn',
'scipy',
'openpyxl',
'tables',
'pillow',
'importlib_metadata'
]

classifiers = [
    'Programming Language :: Python :: 3',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Operating System :: Unix',
    'Operating System :: MacOS :: MacOS X',
    'Intended Audience :: Science/Research'
]

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="plateypus",
    version="0.3.9",
    author = "Sooraj Achar",
    author_email = "acharsr@nih.gov",
    description = "Processes and plots high throughput cytometry experiments through a GUI",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/soorajachar/plateypus",
    classifiers = classifiers,
    packages=find_packages(),
    install_requires=install_requires,
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "plateypus = plateypus.__main__:main"
        ]
    },
)
