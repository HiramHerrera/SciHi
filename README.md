# SciHi

Repository for the analysis of data of SCI-HI observations, these notebooks and  python files were compiled using:

Python - 2.7.15 

##Needed Libraries
| Package | Version |               
|---|---|
| astropy | 2.0.8 |
| h5py | 2.8.0 |
| healpy | 1.12.4 |
| matplotlib | 2.2.3 | 
| numpy | 1.15.1 |
| pandas | 0.23.4 |
| scipy | 1.1.0 |
| seaborn | 0.9.0 |


All of these libraries, except for healpy are within the Anaconda instalation but if needed they may be installed using conda with the usual command:

```
conda install packagename
```

To install healpy just run the following commands:

```
conda config --add channels conda-forge
conda install healpy
```

Also theres a library needed to run these notebooks and codes, which is named [pygsm](https://github.com/telegraphic/PyGSM), if you want to install this library run the following commands:

```
git clone https://github.com/telegraphic/PyGSM \
cd PyGSM
python setup.py install
```

The data needed for these notebooks such as the antenna beam pattern and data will be uploaded further.


