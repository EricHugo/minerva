from setuptools import setup, find_packages

with open('README.rst', mode='r') as f:
    l_description = f.read()

setup(name='minerva',
      version='0.8.0rc1',
      description='Systematic of gene presence/absence and synteny',
      long_description=l_description,
      url='https://github.com/EricHugo/minerva',
      author='Eric Hugoson',
      author_email='eric@hugoson.org',
      license='GNU General Public License v3 or later (G  PLv3+)',
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Environment :: Console',
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          ],
      keywords='bioinformatics comparative genomics',
      install_requires=[
          'biopython>=1.53',
          'numpy>=1.13',
          'matplotlib>=2.0.2',
          'termcolor>=1.1.0',
          'micomplete>=1.1.0',
          ],
      packages=find_packages(),
      include_package_data=True,
      entry_points={
          'console_scripts': [
              'minerva = minerva.minerva:main',
              ]
            }
      )
