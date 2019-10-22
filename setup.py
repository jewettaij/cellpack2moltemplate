from setuptools import setup

setup(
  
  name='cellpack2moltemplate',

  packages=['cellpack2moltemplate'],

  author='Andrew Jewett',

  author_email='jewett.aij@gmail.com',

  url='https://github.com/jewettaij/cellpack2moltemplate',

  download_url='https://github.com/jewettaij/cellpack2moltemplate/archive/v0.1.0.zip',

  version='0.1.0',

  keywords=['CellPACK', 'moltemplate', 'simulation', 'LAMMPS'],

  license='MIT',

  classifiers=['Development Status :: 3 Alpha',
               'Environment :: Console',
               'License :: OSI Approved :: MIT License',
               'Operating System :: MacOS :: MacOS X',
               'Operating System :: POSIX :: Linux',
               'Operating System :: Microsoft :: Windows',
               'Programming Language :: Python',
               'Programming Language :: Unix Shell',
               'Topic :: Multimedia :: Graphics :: 3D Modeling',
               'Topic :: Scientific/Engineering :: Chemistry',
               'Intended Audience :: Science/Research'],

  entry_points={
      'console_scripts': ['cellpack2lt.py=cellpack2moltemplate.cellpack2lt:main']
      },
  zip_safe=True,
  include_package_data=True
)
