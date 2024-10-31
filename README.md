# MATLAB codes for article "Analytical modeling of evanescent coupling in metasurface absorbers for enhanced low-frequency sound control"

These MATLAB codes can be used to calculate the absorbing properties of a planar metasurface absorber composed of detuned resonators arranged in a regular rectangular grid. The resonators are arranged in rectangular super-cells (elementary units) which repeat periodically or with mirror symmetry. The resonators can have rectangular or circular cross-sections and they are characterized by their input specific acoustic impedances.

The main file is **Evanescent.m**. Running this file in MATLAB calculates the absorption coefficients depicted in Figure 4 in the article, namely:
* option **figure = "4(a)"** - rectangular elements, periodic symmetry,
* option **figure = "4(b)"** - rectangular elements, mirror symmetry,
* option **figure = "4(c)"** - circular elements, periodic symmetry,
* option **figure = "4(d)"** - circular elements, mirror symmetry.

The absorption coefficient spectra are compared (depicted together with) the spectra calculated employing the FEM-analytical model.  All the parameters are set in the file
**Evanescent.m** except for the impedances of the individual resonators, which are set in file **Zin.m**.


Author: [Milam Cervenka](https://phys.fel.cvut.cz/en/person/?who=cervenm3&jaz=en), <milan.cervenka@fel.cvut.cz>

