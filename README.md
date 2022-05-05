# Atmospheric Correction
Scripts to modify [SeaDAS](https://seadas.gsfc.nasa.gov/) for using [SS14](https://doi.org/10.1016/j.rse.2013.12.001) and [SSP](https://doi.org/10.1364/OE.27.0A1118) algorithms.

*Note that these scripts are tested with [SeaDAS](https://seadas.gsfc.nasa.gov/) v2020.1 and might need slight modifications to be used with newer versions.*


*To install these scripts in [SeaDAS](https://seadas.gsfc.nasa.gov/), it should be compiled. The process to compile [SeaDAS](https://seadas.gsfc.nasa.gov/) can be found [here](https://seadas.gsfc.nasa.gov/build_ocssw/#building-the-code).*


## How to use these scripts?

* Backup the existing scripts in $OCSSWROOT/ocssw_src/src/l2gen and then place the scripts provided here.
* Add calc_kappa.c to CMakeLists.txt present in $OCSSWROOT/ocssw_src/src/l2gen in both l2gen and l3gen executable sections.
* Compile the code as described [here](https://seadas.gsfc.nasa.gov/build_ocssw/#building-the-code).
* While running l2gen use aer_opt=-20 for [SSP algorithm](https://doi.org/10.1016/j.rse.2013.12.001) and -21 for [SS14 algorithm](https://doi.org/10.1364/OE.27.0A1118).


* If you want to use kappa as a L2 product you will need to modify l2prod.h, l12_proto.h and prodgen.c as mentioned in HOWTO_Add_a_product.txt present in the $OCSSWROOT/ocssw_src/src/l2gen folder inside [ocssw](https://oceancolor.gsfc.nasa.gov/docs/ocssw/index.html). You will have to add the product in the product.xml file present in $OCDATAROOT/common.

## References

* Singh, R.K., Shanmugam, P., 2014. A novel method for estimation of aerosol radiance and its extrapolation in the atmospheric correction of satellite data over optically complex oceanic waters. Remote Sens. Environ. 142, 188–206. https://doi.org/10.1016/j.rse.2013.12.001


* Singh, R.K., Shanmugam, P., He, X., Schroeder, T., 2019. UV-NIR approach with non-zero water-leaving radiance approximation for atmospheric correction of satellite imagery in inland and coastal zones. Opt. Express 27, A1118–A1145. https://doi.org/10.1364/OE.27.0A1118


## Author
**Rakesh Kumar Singh**
