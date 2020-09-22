from setuptools import setup, find_packages


with open('requirements.txt') as f:
    REQUIREMENTS = f.read().strip().split('\n')

with open('VERSION') as f:
    VERSION = f.read().strip()

with open('README.rst') as f:
    README = f.read()

setup(
    name='rootwater',
    version=VERSION,
    license='MIT',
    description='Python package to estimate root water uptake (RWU) from rhizosphere soil moisture dynamics and some tools to handle sap flow measurements',
    long_description=README,
    url="https://github.com/cojacoo/rootwater",
    author='Conrad Jackisch',
    author_email='conrad.jackisch@tbt.tu-freiberg.de',
    install_requires=REQUIREMENTS,
    test_require=['nose'],
    test_suite='nose.collector',
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    zip_safe=False
)
