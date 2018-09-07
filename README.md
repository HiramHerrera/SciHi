# SciHi

Repository for the analysis of data of SCI-HI observations, these notebooks and  python files were compiled using:

[Python](https://www.python.org/) - 2.7.15 

## Needed Libraries

| Package | Version |               
|---------|:--------|
| [astropy](http://www.astropy.org/) | 2.0.8 |
| [h5py](https://www.h5py.org/) | 2.8.0 |
| [healpy](https://healpy.readthedocs.io/en/latest/) | 1.12.4 |
| [matplotlib](https://matplotlib.org/) | 2.2.3 | 
| [numpy](http://www.numpy.org/) | 1.15.1 |
| [pandas](https://pandas.pydata.org/) | 0.23.4 |
| [scipy](https://www.scipy.org/) | 1.1.0 |
| [seaborn](https://seaborn.pydata.org/) | 0.9.0 |


All of these libraries, except for healpy are within the Anaconda instalation but if needed they may be installed using conda with the usual command:

```
conda install packagename
```

To install healpy just run the following commands:

```
conda config --add channels conda-forge
conda install healpy
```

Also to run these codes and notebooks you will need [pygsm](https://github.com/telegraphic/PyGSM), if you want to install this library run the following commands:

```
git clone https://github.com/telegraphic/PyGSM
cd PyGSM
python setup.py install
```

The data needed for these notebooks such as the antenna beam pattern and data will be uploaded further.


