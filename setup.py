from setuptools import setup

def readme():
      with open('README.rst') as f:
               return f.read()

setup(name='TEsingle',
      version='1.0',
      description='Tool for estimating differential enrichment of Transposable Elements and other highly repetitive regions in single-cell data',
      long_description_content_type="text/x-rst",
      long_description=readme(),
      classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'Environment :: Console',
            'Natural Language :: English',
            'License :: OSI Approved :: BSD License',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Programming Language :: Python :: 3',
            'Operating System :: Unix'
      ],
      keywords='TE transposable element differential enrichment single-cell',
      url='https://www.mghlab.org/software/tesingle',
      author='Talitha Forcier, Cole Wunderlich, Oliver Tam, Molly Gale Hammell',
      author_email='mghcompbio@gmail.com',
      license='BSD-3-Clause',
      packages=[
          'TEtools'
      ],
      platforms=[
          'Linux'
      ],
      install_requires=[
          'networkx',
          'scipy',
          'numpy',
          'pysam>=0.9'
      ],
      include_package_data=True,
      zip_safe=False,
      scripts=[
          'TEsingle'
      ]
)
