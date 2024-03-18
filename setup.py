from setuptools import setup

def readme():
      with open('README.rst') as f:
               return f.read()

setup(name='TEsingle',
      version='0.1',
      description='Tool for estimating differential enrichment of Transposable Elements and other highly repetitive regions in single-cell data',
      long_description=readme(),
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Environment :: Console',
          'Natural Language :: English',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Programming Language :: Python :: 3',
          'Operating System :: Unix'
      ],
      keywords='TE transposable element differential enrichment single-cell',
      url='http://hammelllab.labsites.cshl.edu/software',
      author='Talitha Forcier, Cole Wunderlich, Oliver Tam, Molly Gale Hammell',
      author_email='talitha.forcier@nyulangone.edu',
      license='GPLv3',
      packages=[
          'TEtools'
      ],
      platforms=[
          'Linux'
      ],
      install_requires=[
          'argparse',
          'pysam>=0.9'
      ],
      include_package_data=True,
      zip_safe=False,
      scripts=[
          'TEsingle'
      ]
)
