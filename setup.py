from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


# Get the program version from another file.
exec(open('porechop/version.py').read())


setup(name='Kleborate',
      version=__version__,
      description='Kleborate',
      long_description=readme(),
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GPLv3',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Intended Audience :: Science/Research',
      ],
      keywords='microbial genomics sequence typing',
      url='https://github.com/katholt/Kleborate',
      author='Kathryn Holt',
      author_email='',
      license='GPLv2',
      packages=['kleborate'],
      install_requires=[
      ],
      test_suite='nose.collector',
      tests_require=[],
      entry_points={
          'console_scripts': ['kleborate = kleborate.kleborate:main'],
      },
      include_package_data=True,
      zip_safe=False)
