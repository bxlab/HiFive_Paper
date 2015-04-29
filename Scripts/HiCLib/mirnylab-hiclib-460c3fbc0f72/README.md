hiclib
======

Read the [documentation](http://mirnylab.bitbucket.org/hiclib/index.html).

Installation
------------
The easiest way is to use pip, but see notes below:

First, install [mirnylib](https://bitbucket.org/mirnylab/mirnylib):

`$ pip install https://bitbucket.org/mirnylab/mirnylib/get/tip.tar.gz`

Then install hiclib:

`$ pip install https://bitbucket.org/mirnylab/hiclib/get/tip.tar.gz`

Alternatively, there is a linux install script that doesn't use pip.

Installation requirements
-------------------------

### Non-python dependencies

- Mirnylib provides a HDF5 dictionary type based on the Python h5py package. For the h5py package to install properly, you need to have the shared library and development headers for HDF5 1.8.4 or newer installed (`libhdf5-dev` or similar). See the [h5py docs](http://docs.h5py.org/en/latest/build.html) for more information.

- You will also need the `bowtie2` mapping software which can be downloaded manually or installed via the download_bowtie.sh script in the hiclib folder.

### Python dependencies

Required:

Getting the basic Scientific Python stack (numpy/scipy/matplotlib) can be trickier on some platforms than others. For more details, see the [instructions on scipy.org](http://www.scipy.org/install.html). You should already have these dependencies installed and running correctly before attempting to install this package.

- numpy (1.6+)
- scipy
- matplotlib

Installed by setuptools if missing:

- biopython
- pysam

Optional:

- Sphinx (to generate documentation)