# To do list

[ ] Improve AR1 estimation
[ ] Boundary rescaling for all kernels (only gaussian currenlty)
[ ] Disable profiling
[ ] Instead of rescaling, change gradient so it is scaled by the total weight
  - make sure hessian also is scaled by the same thing


# Version 0.0.2

* Updated computation of EBIC (was previously wrong)
* Updated computation of profile likelihood in covariance update (previously was not using dispersion)
* Made Kernel abstract and added Epanechnikov kernel
* Made WorkingCovariance abstract and Independent and added AR(1)
* AR(1) estimation is not currently working.
* Improve stepsize adaptation and dropped backtracking line search
* Added sine function in data generation
* Data generation now allows user-supplied group difference function
* Made penalty abstract and added SCAD and MCP
* Uniformize initialization to independent MLE
* Added rescaling of the kernel at the boundary
* Initialize covariance parameters under unpenalized independent model
* Enable two-step estimation where covariance parameters are fixed after initialization
* Added profiling to identify expensive functions

# Version 0.0.1

* First release
* Added documentation and resolved R CMD CHECK warnings and notes.

# Version 0.0.0.9002

* Version for JSM poster

# Version 0.0.0.9001

* C++ backend implemented
* Data generating example implemented

# Version 0.0.0.9000

* Package initialization.
