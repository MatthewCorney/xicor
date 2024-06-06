from setuptools import setup, find_packages

# Load the README file as the long description
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='xicor',
    version='0.0.1',
    description="Numerically equivalent reimplementation of the the R xicor package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/MatthewCorney/xicor',
    author='Matthew Corney',
    author_email='matthew_corney@yahoo.co.uk',
    package_dir={'': "src"},
    packages=find_packages(where="src"),
    install_requires=[
        'numpy', 'scipy',
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    extras_require={
        'dev': ['pytest', 'pytest-runner'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    zip_safe=False
)