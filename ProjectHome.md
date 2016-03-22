From Monty Python and the Holy Grail
  * LAUNCELOT:  Look, my liege!
  * ARTHUR:  Camelot!
  * GALAHAD:  Camelot!
  * LAUNCELOT:  Camelot!
  * PATSY:  It's only a model.
# News #
  * **16 July 2010**: CAMELOT GUI, version 0.2 is released. The new features include scrollbar for the tabs, additional utilities to compute brine density and viscosity, and a tab for Area of Review calculations.
  * **10 May 2010**: CAMELOT GUI is released
  * **13 January 2010**: initial release of numerical engines

# Overview #
CAMELOT is a collection of Python modules to model the injection induced brine pressure increase. The program is intended to model injection of carbon dioxide into deep saline aquifers as part of geologic carbon sequestration using simple semi-analytic solutions. The modules can either be directly integrated into scripts, or they can be run from the CAMELOT GUI.

## Numerical engines ##
Currently three numerical modules for injection into a single aquifer are implemented. For all three it is assumed that the aquifer is homogeneous and isotropic, and that an equivalent volume of brine is injected instead of CO<sub>2</sub> (to allow for single-phase modeling). The three modules are:
  * **[Theis](Theis.md)**: injection into an aquifer with infinite lateral extent bound at the top and bottom by impermeable layers. . The well radius is assumed to be infinitesimally small. Requires the module [SciPy](http://www.scipy.org/).
  * **Hantush-Jacob**: injection into an aquifer with infinite lateral extent bound at the top and bottom by leaky layers. The well radius is assumed to be infinitesimally small, and storage in the leaky layers is neglected. Requires the module [SciPy](http://www.scipy.org/).
  * **Moench with Zhou et al. extension** (SingleLayer): injection into an aquifer either of infinite lateral extent or with a lateral impermeable boundary bound at the top and bottom by leaky layers. In contrast to the Hantush-Jacob solution the well may have a finite radius and storage in the leaky layers is not neglected. Requires the module _HoogPy_ for numerical inversion of the Laplace transform and the modules [SciPy](http://www.scipy.org/) and [NumPy](http://numpy.scipy.org/).
## Utilities ##
CAMELOT comes with several utility modules:
  * **Utilities**: includes functions to compute hydraulic conductivity, hydraulic transmissivity, specific storage, storativity, volumetric injection rate of equivalent brine, and initial aquifer pressure.
  * **Units**: enables unit conversion.
  * **CO2Properties**: computes the density and dynamic viscosity of CO<sub>2</sub> based on pressure and temperature.
  * **BrineProperties**: computes the density and dynamic viscosity of brine based on temperature, pressure and salinity.
## GUI ##
All the above modules can be accessed through the CAMELOT GUI. The GUI allows to select both aquifer properties and injection parameters; runs the numerical engines; plots the results; and exports the results in csv format. The GUI is invoked by executing the _runGUI_ script. The module _CAMELOT_ contains the GUI, and additional widgets used in the GUI are stored in the module _ExtraWidgets_.
# Installation #
The GUI and the necessary numerical engines and utilities are installed either using the executable installer (Windows only) or by using [distutils](http://docs.python.org/distutils/). The third-party packages [NumPy](http://numpy.scipy.org/), [SciPy](http://www.scipy.org/), and [matplotlib](http://matplotlib.sourceforge.net/) need to be installed separately.
  * To use the Windows installer, download the executable to any directory and run the installer by double clicking it. This will copy all necessary modules to _Lib/site-packages/CAMELOT_ in your Python directory, as well as creating the _CAMELOT.pth_ file. The script that starts the GUI, _runCAMELOT\_X\_X.py_ needs to be downloaded separately, and can be run from any folder.
  * To run [distutils](http://docs.python.org/distutils/) directley the zip file _CAMELOT-X.X.zip_ needs to be downloaded and unpacked into any folder. The [distutils](http://docs.python.org/distutils/) command _install_ is then run from the command line: `>python setup.py install`. This will copy all necessary modules to _Lib/site-packages/CAMELOT_ in your Python directory, as well as creating the _CAMELOT.pth_ file. The script that starts the GUI, _runCAMELOT\_X\_X.py_ needs to be downloaded separately, and can be run from any folder.
# Downloads #
The CAMELOT GUI can be downloaded using an installer. The modules that are used in the GUI can also be downloaded separately.