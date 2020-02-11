from setuptools import setup, find_packages
setup(
	name='ATL11',
	version='1.0.0.1',
	description='Produces the ICESat-2 Land-Ice Along-Track H(t) product (ATL11)',
	url='https://github.com/suzanne64/ATL11',
	author='Ben Smith',
	author_email='besmith@uw.edu',
	license='MIT',
	classifiers=[
		'Development Status :: 3 - Alpha',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Physics',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 2',
		'Programming Language :: Python :: 2.7',
	],
	keywords='altimetry, spatial queries, ICESat-2',
	packages=find_packages(),
        include_package_data=True,
	install_requires=['numpy','matplotlib','h5py','gdal','scipy'],
)
