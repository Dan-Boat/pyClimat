from setuptools import setup
import sys
import pathlib

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setup(
      name = "pyClimat",
      version = "0.1.0",
      description= "Python package for climate model output analysis and visualization",
      long_description = README,
      long_description_content_type = "test/markdown",
      keywords = "climate analysis, climate data visualization, eof, Cartopy, Xarray",
      url="https://github.com/Dan-Boat/PyClimat",
      author="Daniel Boateng",
      author_email= "dannboateng@gmail.com",
      license="MIT",
      packages=["pyClimat"],
      install_requires=["mpi4py",
                        "numpy", 
                        "pandas",
                        "xarray", 
                        "statsmodels",
                        "seaborn",
                        "sklearn",
                        "scipy",
                        "eofs",
                        "Cartopy",
                        "matplotlib" 
                        ],

      classifiers = [
        'Development Status :: 1 -',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5', 
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8', 
        'Programming Language :: Python :: 3.9',
        
          ],
      
      include_package_data = True,
      entry_points = {}
      
      )