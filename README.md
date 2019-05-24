# DE-LA-MO

The De-la-mo website can be found at http://thermal.cnde.iastate.edu/de-la-mo.xhtml

On-line documentation and installation instructions are located at
http://thermal.cnde.iastate.edu/de-la-mo/index.html

De-La-Mo is Copyright (C) 2016-2018 Iowa State University
Research Foundation, Inc. It is published under the
Apache 2.0 license. See the LICENSE.txt file for details


## Repository Organization
* `doc/': Contains documentation source code
* `scripts/': Contains Python helper scripts for building De-la-mo models
* `delamo/`: Contains the python package _delamo_
* `examples/`: Contains example python scripts which demonstrate how to use _delamo_
* `examples/data/`: Contains data files used in the example scripts
* `autofiberlib/`: Git submodule to the [AutoFiberlib](https://github.com/nscheirer/autofiberlib) repository. Contains
advanced functions for optimizing composite fiber layups based on geodesic paths and strain energy minimization. (Optional)

## Quickstart Build Instructions
1. Acquire requirements listed at http://thermal.cnde.iastate.edu/de-la-mo.xhtml
2. python setup.py build
3. python setup.py install

## Cloning Instructions
In order to include the autofiber capabilities during cloning the following command should be run:
```
git clone --recursive [Url to this repo]
```

If you have already cloned this repository the following is used to load the autofiberlib submodule:
```
git submodule update --init
```

## Installation Structure

## How to Run De-la-mo

Each of the Python scripts in the `examples/' directory can be run directly
from a Python interpreter. This will create a subdirectory with the
`_output' suffix containing the generated CAD model (`.step' file) and
generated Python script for ABAQUS (`.py' file). Run the generated `.py' file
inside ABAQUS/CAE to perform the desired analysis. 

The 'delamo_process' utility (in the `scripts/' directory) can be used as a
helper to assist in running De-la-mo scripts. It is especially helpful
for multi-step processes such as automated damage insertion.
