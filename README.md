# The Vaganov-Shashkin tree-ring growth model (VSM) - version Oscilloscope

The repository is based on the version of the Vaganov-Shashkin model published by Anchukaitis et al. (2020) available through http:/github.com/kanchukaitis/vsm/. The referenced MATLAB code fully reproduced the original FORTRAN implementation of the Vaganov-Shashkin model. However, because modified versions (mostly simplifications) of this overal complex tree-ring model were published over time, it does not produce results consistent with all currently existing model implementations. Here we introduce modifications making the calculus consistent with the version of the model used in the visual tool VS-Oscilloscope v1.367 (http://www.vs-genn.ru/).

The modifications are based on introduction of two new functions (`vsm_oscilloscope` and `NEXT`), we also provide new file showing the required structure of parameters input (`parameters_oscillosope`). The original functions were not modified at all. With the set of these tree new files, the modified version of the model is self-standing and might be used independently from functions used in original repository (although we want to fully account for the fact that the development of the new functions was largely permited by the existence of original function `vsm` and script `parameters` which we modified). The details about application syntax is provided directly in the heading of `vsm_oscilloscope` function, where each modification made to the original `vsm` function was also annotated.

## Citation

To refer `vsm` and additional functions from the original repository (which served as baseline for development of Oscillosope versions) please cite:

Anchukaitis K.J., Evans M.N., Hughes M. K.,  Vaganov E.A. (2020): An interpreted language implementation of the Vaganov-Shashkin tree-ring proxy system model. *Dendrochronologia* 60: 125677. DOI 10.1016/j.dendro.2020.125677

More details about this implementation are available through http:/github.com/kanchukaitis/vsm/

The modified function `vsm_oscilloscope` was recently presented in following publications:
