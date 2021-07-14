# The Vaganov-Shashkin tree-ring growth model (VSM) - version *Oscilloscope*

The repository is based on the version of the Vaganov-Shashkin model published by Anchukaitis et al. (2020) available through https://github.com/kanchukaitis/vsm/. The referenced MATLAB code fully reproduces the original FORTRAN implementation of the Vaganov-Shashkin model. However, because modified versions (mostly simplifications) of this overall complex tree-ring model were published over time, it does not provide consistent results with all currently existing model implementations. Here we introduce modifications mostly to the original `vsm` function making the calculus consistent with the version of the model used in the visual tool VS-Oscilloscope v1.367 (http://www.vs-genn.ru/).

The modifications are based on the introduction of two new functions (`vsm_oscilloscope` and `NEXT`), we also provide a new file showing the required structure of input parameters (`parameters_oscilloscope`). The original functions were not modified at all and, thus, can be applied fully in the manner described by Anchukaitis et al. (2020). The modified version of the model is self-standing and might be used independently from functions used in the original repository (although we want to fully account for the fact that the syntax of the new functions was strongly inspired by original functions `vsm` and script `parameters`). The details about application syntax are provided directly in the heading of `vsm_oscilloscope` function, where each modification made to the original `vsm` function was also annotated.

## Citation

To refer `vsm` and additional functions from the original repository (which served as baseline for the development of Oscilloscope versions) please cite:

   - Anchukaitis K.J., Evans M.N., Hughes M. K.,  Vaganov E.A. (2020): An interpreted language implementation of the Vaganov-Shashkin tree-ring proxy system model. *Dendrochronologia* 60: 125677. https://doi.org/10.1016/j.dendro.2020.125677
   
   - For additional details and up-to-date version of this implementation please refer to https://github.com/kanchukaitis/vsm/

The modified function `vsm_oscilloscope` was recently presented in the following publications:

   - Tumajer J., Kaspar J., Kuzelova H., Shishov V.V., Tychkov I.I., Popkova M.I., Vaganov E.A., Treml V. (2021): Forward Modeling Reveals Multidecadal Trends in Cambial Kinetics and Phenology at Treeline. *Frontiers in Plant Science* 613646. https://www.frontiersin.org/articles/10.3389/fpls.2021.613643/full
   - Tumajer J., Buras A., Camarero J.J., Carrer M., Shetti R., Wilmking M., Altman J., Sangüesa-Barreda G., Lehejcek J.: Growing faster, longer or both? Modelling plastic response of *Juniperus communis* growth phenology to climate change. Submitted to *Global Ecology and Biogeography.*
   - Tumajer J., Shishov V.V., Ilyin V.A., Camarero J.J.: Plastic growth dynamics of Mediterranean pines and junipers determines their climatic adaptability. Submitted to *Agricultural and Forest Meteorology.*

## User notes

The newly-introduced Oscilloscope functions were fully developed in Octave and were NOT tested in MATLAB. However, the code is build exclusivelly using core Octave functions and, thus, proper functionality with MATLAB might be assumed. Please, do not hesitate to contact me if you spot any problems.

## Acknowledgements

The identification of the differences between both VS-model implementations and the development of modified functions was performed as a part of following projects:

   - *Alexander von Humboldt Foundation - PostDoctoral Fellowship to Jan Tumajer*
   - *Czech Science Foundation - project 19-138076S*
   - *Charles University Center for Excellence - project UNCE/HUM 018*

## Contact and bug-reporting
tumajerj@natur.cuni.cz
