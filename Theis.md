#description of Theis module.

# Introduction #

The _Theis_ module implements the solution published by Theis (1935) to the problem of injection into an infinite aquifer bounded by impermeable layers at the top and bottom. The solution assumes an infinitesimal well radius, and a homogenous and isotropic aquifer.

# Model set up #
Three parameters are needed to set up the model:
  * aquifer transmissivity (T)
  * aquifer storativity (S)
  * volumetric injection rate (Q)
These three parameters are set using the function `SetRequiredProperties(Q, S, T)`. In addition several additional parameters can be specified by adding optional parameters to the function call:
  * InjectionEnd: the time at which injection stops (default: continuous injection)
  * HeadOrPressure: report result as head (H) or pressure (P) (default: head)
  * FresWDensity: density of fresh water; used to translate head to pressure (default: 1000)
  * Gravity: gravitational acceleration (default: 9.81)
  * WarnNoRadius: warning message that is returned if no radius could be computed for the given time and head/pressure (default: the string "no contour")
  * WarmGTsteadystate: warning message that is returned if no time could be computed for the given radius and head/pressure (default: the string "pressure not reached")
All these parameters (required and optional) can also be set by the function `SetProperty(propertytype, value)`.

# Model calculations #
A Theis model has three variables: the radius, the time, and the head/pressure. If two of these variables are given, the third can be calculated. Thus there are three functions to compute results:
  * **GetHeadForRadiusTime(radius, time)**: computes either head or pressure based on the given radius and time, and on the required properties. Radius and time can be either single values or lists.
  * **GetRadiusForHeadTime(head/pressure, time)**: computes the radius based on the given head or pressure and time, and on the required properties. Head/pressure and time can be either single values or lists. If no radius can be found the warning message specified by _WarnNoRadius_ is returned.
  * **GetTimeForRadiusHead(radius, head/pressure)**: computes the time based on the given head or pressure and radius, and on the required properties. Head/pressure and radius can be either single values or lists. If no time can be found the warning message specified by _WarnGTsteadystate_ is returned.

# Units #
The _Theis_ module can be used with any set of consistent units (e.g., all length meassures in meters, all time meassures in seconds, ...).

# Documentation #
For more detailed instructions on how the _Theis_ module is used please refer to the documentation on the download page.

# Example #
For an example on how to use the _Theis_ module please refer to the _ExampleTheis\_X\_X.py_ script on the download page.

# References #
  * C. V. Theis (1935). The Relation Between Lowering of the Piezometric Surface and the Rate and Duration of Discharge of a Well Using Ground Water Storage. Transactions American Geophysical Union, 16th anual meeting, Pt. 2, pp. 519-524