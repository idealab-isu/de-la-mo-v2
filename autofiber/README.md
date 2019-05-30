# AutoFiber
Strain energy minimization of geodesic based parameterizations over 3D
surfaces for optimization of fiber layup orientations. Created to be part
of the [De-La-Mo](http://thermal.cnde.iastate.edu/de-la-mo.xhtml)
automatic defect insertion into FEM package developed at Iowa
State University.

AutoFiber is Copyright (C) 2016-2018 Iowa State University
Research Foundation, Inc. It is published under the
Apache 2.0 license. See the LICENSE file for details.

## Project Structure
* `autofiber/`: Contains the python package __autofiber__.
* `demos/`: Contains a variety of demo models and a script, *test.py*
which demonstrates usage on each model.
* `doc/`: Contains various documentation materials

## Package Structure
* `generator.py`: Main control script which generates geodesic start points, calculates geodesic paths, computes the
geodesic parameterization, and employs the strain energy minimization.
* `geodesic.py`: Contains the functions necessary for computing the geodesic trajectories and parameterizations.
* `optimization.py`: Contains the definition for the strain energy function and gradient as well as the RMSprop
algorithm used to minimize the strain energy function.
* `analyze_uv.py`: Spatialnde helper script that defines mesh adjacency indexes.

## Dependencies
* Tested with Python 2.7-3.7
* `spatialnde`: 3D model loader and image projection package \
Created by Dr. Stephen D. Holland at Iowa State University \
[Spatialnde](http://thermal.cnde.iastate.edu/spatialnde)
* `Numpy`
* `Matplotlib` - optional for plotting

## Installation
Once all dependencies are installed run:
```
python setup.py build
python setup.py install
```

## How to run
Take a look at `test.py` for an in-depth explanation of the relevant
API calls and how they work for a variety of models.

## Abaqus integration
A function has been created and implemented in De-La-Mo to allow for
automatic insertion of computed fiber orientations into Abaqus. See the
[De-La-Mo](http://thermal.cnde.iastate.edu/de-la-mo.xhtml) project for
more information.
