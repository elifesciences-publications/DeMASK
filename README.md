# DeMASK
## **De**terminstic **M**odeling for **A**nalysis of complex **S**ingle molecule **K**inetics

This package was designed to assist those attempting to understand complicated single molecule kinetics by arming them with tools to test kinetic models and fit microscopic rate constants. To use this tool, you must perform dwell time analysis and generate cumulative distributions of dwell times.

* MATLAB software for deterministic modeling and fitting with single molecule data 

## Related Resources
* <a href="https://www.biorxiv.org/content/early/2018/05/10/319749">Our preprint</a>

## Included Content
* DeMASK.m - main MATLAB fitting/bootstrap script
* smODEfcnON.m - function defining ODEs
* myresODEgloOnOffFull.m - function generating residuals used in least squares minimization
* myresODEgloFIGSon.m - function to generate figures after fit

## Instructions
* Download repository and store all contects together in one folder (all .m MATLAB scripts and all .csv data files).
* Run DeMASK script. 
* When prompted, enter the path to the scripts/data.
* When prompted enter the number of desired bootstrap iterations to perform.
