import pathlib
from setuptools import setup
from Cython.Build import cythonize

setup(
    name='msmutect',
    version='4.1',
    packages=['src', 'src.Entry', 'src.GenomicUtils', 'src.IndelCalling'],
    entry_points={
        'console_scripts': ['msmutect=src.Entry.main:main'],
    },
    url='https://github.com/MaruvkaLab/MSMuTect_0.5',
    license='MIT',
    author='Avraham Kahan, Yossi Maruvka, and the Maruvka Lab at Technion',
    author_email='yosi.maruvka@bfe.technion.ac.il',
    description='Tool to determine microsatellite in/stability in Tumors from DNA sequencing',
    long_description=f"{pathlib.Path(__file__).parent}/README.md",
    long_description_content_type="text/markdown",
    install_requires=['typing>=3.7.4.3',
                        'numpy>=1.20.1',
                        'pysam>=0.16.0.1',
                        'scipy>=1.6.1',
                        'setuptools>=54.1.1'
                        ],

    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    ext_modules=cythonize([
                        'src/IndelCalling/FisherTest.pyx',
                        'src/IndelCalling/CallMutations.pyx',
                         'src/IndelCalling/CallAlleles.pyx',
                         'src/IndelCalling/Histogram.pyx',
                        'src/IndelCalling/AlleleSet.pyx',
                         'src/IndelCalling/Locus.pyx',

                        'src/GenomicUtils/ReadsFetcher.pyx',
                         'src/GenomicUtils/LocusFile.pyx',

                           'src/Entry/PairFileBatches.pyx',
                           'src/Entry/SingleFileBatches.pyx'],
                          compiler_directives={'language_level' : "3", "profile": True})
)
