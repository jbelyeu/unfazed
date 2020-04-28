import re

from setuptools import find_packages, setup


with open("README.md", "r") as fh:
    long_description = fh.read()

with open("unfazed/__init__.py", "r") as fd:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', fd.read(), re.MULTILINE
    ).group(1)

with open("requirements.txt", "r") as f:
    requires = f.read().splitlines()


setup(
    name="unfazed",
    version=version,
    description="command line tool for extended read-backed phasing of SNVs and SVs, plus allele-balance phasing for CNVs",
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Jonathan Belyeu",
    author_email="jrbelyeu@gmail.com",
    url="https://github.com/jbelyeu/unfazed.git",
    packages=["unfazed"],
    package_data={"": ["LICENSE", "README.md"]},
    include_package_data=True,
    install_requires=requires,
    license="MIT",
    zip_safe=False,
    entry_points={"console_scripts": ["unfazed = unfazed.__main__:main"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
