#!/usr/bin/env python

from setuptools import setup, find_packages
import versioneer

setup(name='graynet',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Individualized single-subject networks from T1 mri features '
                  'such as cortical thickness and gray matter density. ',
      long_description="""
Individualized single-subject networks from T1-weighted MRI features such as:
  - Cortical thickness
  - Gray matter density
  - Subcortical morphometric features
  - Gyrification and curvature

Applicable for whenever network-level features are useful, among which common use cases are: 
 - Biomarker development.
 - Brain-behaviour relationships (e.g. for the diagnosis and prognosis of many brain disorders such as Alzheimer's, Parkinson's, Schizophrenia and the like).
 - Aging (changes in network properties over age and their relations to other variables).

Docs: https://raamana.github.io/graynet/
""",
      author='Pradeep Reddy Raamana',
      author_email='raamana@gmail.com',
      url='https://github.com/raamana/graynet',
      packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
      install_requires=['numpy', 'hiwenet', 'pyradigm', 'nibabel', 'networkx',
                        'medpy', 'scipy'],
      classifiers=[
              'Intended Audience :: Science/Research',
              'Programming Language :: Python',
              'Topic :: Software Development',
              'Topic :: Scientific/Engineering',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS',
              'Programming Language :: Python :: 3.6',
          ],
      entry_points={
          "console_scripts": [
              "graynet=graynet.__main__:main",
          ]
      },
      package_dir ={'graynet':'graynet'},
      package_data={'graynet':['atlases/*', 'resources/*']},
      include_package_data=True

     )
