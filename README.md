# The Vaganov-Shashkin tree-ring growth model (VSM) - version Oscilloscope

The repository is based on the version of the Vaganov-Shashkin model published by Anchukaitis et al. (2020) available through https://github.com/kanchukaitis/vsm/. The referenced MATLAB code fully reproduced the original FORTRAN implementation of the Vaganov-Shashkin model. However, because modified versions (mostly simplifications) of this overal complex tree-ring model were published over time, it does not produce results consistent with all currently existing model implementations. Here we introduce modifications making the calculus consistent with the version of the model used in the visual tool VS-Oscilloscope v1.367 (http://www.vs-genn.ru/).

The modifications are based on introduction of two new functions (`vsm_oscilloscope` and `NEXT`), we also provide new file showing the required structure of parameters input (`parameters_oscillosope`). The original functions were not modified at all. With the set of these tree new files, the modified version of the model is self-standing and might be used independently from functions used in original repository (although we want to fully account for the fact that the development of the new functions was largely permited by the existence of original function `vsm` and script `parameters` which we modified). The details about application syntax is provided directly in the heading of `vsm_oscilloscope` function, where each modification made to the original `vsm` function was also annotated.

## Citation

1. To refer `vsm` and additional functions from the original repository (which served as baseline for development of Oscillosope versions) please cite:

   - Anchukaitis K.J., Evans M.N., Hughes M. K.,  Vaganov E.A. (2020): An interpreted language implementation of the Vaganov-Shashkin tree-ring proxy system model. *Dendrochronologia* 60: 125677. https://doi.org/10.1016/j.dendro.2020.125677
   
   - More details about this implementation are available through https://github.com/kanchukaitis/vsm/

The modified function `vsm_oscilloscope` was recently presented in following publications:

   - Tumajer J., Kaspar J., Kuzelova H., Shishov V.V., Tychkov I.I., Popkova M.I., Vaganov E.A., Treml V. (2021): Forward Modeling Reveals Multidecadal Trends in Cambial Kinetics and Phenology at Treeline. *Frontiers in Plant Science* 613646. https://www.frontiersin.org/articles/10.3389/fpls.2021.613643/full
   - Tumajer J., Buras A., Camarero J.J., Carrer M., Shetti R., Wilmking M., Altman J., Sanguessa-Bareda G., Lehejcek J.,: Growing faster, longer or both? Modelling plastic response of Juniperus communis growth phenology to climate change. Submitted to *Global Ecology and Biogeography.*
   - Tumajer J., Shishov V.V., Ilyin V.A., Camarero J.J.: Plastic growth dynamics of Mediterranean pines and junipers determines their climatic adaptability . Submitted to *Agricultural and Forest Meteorology.*

