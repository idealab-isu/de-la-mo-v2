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
* `autofiber/`: The [AutoFiber](https://github.com/nscheirer/autofiber) python package. Which contains a module of
advanced functions for optimizing composite fiber layups based on geodesic paths and strain energy minimization.

## Quickstart Build Instructions
1. Acquire requirements listed at http://thermal.cnde.iastate.edu/de-la-mo.xhtml
2. python setup.py build
3. python setup.py install

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

## Updating AutoFiber
Autofiber is a subpackage of De-la-mo, therefore we have implemented a subtree system for including the Autofiber
package. See [here](https://git-scm.com/book/en/v1/Git-Tools-Subtree-Merging) for more details. In order to update
Autofiber run the following commands (make sure you are in the root directory of de-la-mo):

If you have never updated Autofiber before you must add the relevant repository and branch as follows:
```
git remote add autofiber_remote https://github.com/nscheirer/autofiber.git
git fetch autofiber_remote
git checkout -b autofiber_branch autofiber_remote/master
```

Now you can switch between `master` (which contains de-la-mo) and `autofiber_branch` (which contains the autofiber
package). We will need to tell git which directory autofiber is contained within like so:

```
git checkout master
git read-tree --prefix=autofiber/ -u autofiber_branch
```

If git complains about the `read-tree` command then the directory may already be linked, especially if autofiber already exists, so you can continue.

Now, if there is an update to autofiber we can switch to the autofiber branch and pull from autofiber_remote. Initially,
you may need to commit the new autofiber package before you can merge or else the package might appear in a weird place.

```
git checkout autofiber_branch
git pull
```

Finally, we need to merge the autofiber branch with the de-la-mo branch as a subtree as follows:

```
git checkout master
git merge --squash --allow-unrelated-histories -s subtree --no-commit autofiber_branch
```

The git merge command will pull in the autofiber branch while not preserving the history of autofiber and squashing any
commits.
