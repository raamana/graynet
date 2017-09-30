#!/usr/bin/env python

from setuptools import setup, find_packages
import versioneer

setup(name='graynet',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Individualized single-subject networks from T1 mri features such as cortical thickness and gray matter density. ',
      long_description='Individualized single-subject networks from T1 mri features such as cortical thickness and gray matter density; graynet',
      author='Pradeep Reddy Raamana',
      author_email='raamana@gmail.com',
      url='https://github.com/raamana/graynet',
      packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]), # ['neuropredict'],
      install_requires=['numpy', 'pyradigm', 'nibabel', 'networkx', 'medpy'],
      classifiers=[
              'Intended Audience :: Science/Research',
              'Programming Language :: Python',
              'Topic :: Software Development',
              'Topic :: Scientific/Engineering',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS',
              'Programming Language :: Python :: 2.7',
          ],
      entry_points={
          "console_scripts": [
              "graynet=graynet.__main__:main",
          ]
      },
      package_dir ={'graynet':'graynet'},
      package_data={'graynet':'atlases/*'},
      include_package_data=True

     )
