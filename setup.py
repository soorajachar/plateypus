#! /usr/bin/env python3
from setuptools import setup,find_packages
    
install_requires = [
'pandas',
'seaborn',
'scipy',
]

classifiers = [
    'Programming Language :: Python :: 3',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Operating System :: Unix',
    'Operating System :: MacOS :: MacOS X',
    'Intended Audience :: Science/Research'
]

setup(
    name="plateypus",
    version="0.1.6",
    author = "Sooraj Achar",
    author_email = "acharsr@nih.gov",
    description = "Processes and plots high throughput cytometry experiments through a GUI",
    long_description = "file: README.md",
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
