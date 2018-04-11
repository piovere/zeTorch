# zeTorch
Analysis of plasma torch data

We produce spectra from feeding various samples into a plasma torch. The samples are mixed in real time and fed into the sample chamber. **DETECTOR DATA TO FOLLOW**. The goal of this project is to produce a neural network classifier to identify the components of the plasma feed material from the input spectra.

# ToDo
- [ ] Make data files that match our spectra containing which Analytes were present
- [ ] Develop input processing for our spectra files (slice out metadata, determine which columns are needed)
- [ ] Pair input and output data sets
- [ ] Input processing and normalization (PCA?)

# Useful Links
- [LMFIT](https://lmfit.github.io/lmfit-py/): Advanced arbitrary curve fitting library. Wrapper around scipy.optimize
- [Akaike information criterion (AIC)](https://en.wikipedia.org/wiki/Akaike_information_criterion): scoring for a goodness of model fit
- [Bayesian information criterion](https://en.wikipedia.org/wiki/Bayesian_information_criterion): Another tool for selecting between competing models
- [Voigt profile](https://en.wikipedia.org/wiki/Voigt_profile): The curve we will be fitting to for the peaks we detect. Combines doppler broadening and Lorentzian interference. Python implementation [here](https://scipython.com/book/chapter-8-scipy/examples/the-voigt-profile/)
- [Example of 1-D signal smoothing in numpy](http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html)
- [Numerical Methods in Fortran](https://websites.pmc.ucsc.edu/~fnimmo/eart290c_17/NumericalRecipesinF77.pdf): Textbook on implementing various numerical methods
- [Voigt profile function in Python](https://scipython.com/book/chapter-8-scipy/examples/the-voigt-profile/)
- [Ocean Optics spectrometer library for Python](https://github.com/ap--/python-seabreeze)
- [Convolution example](https://stackoverflow.com/questions/40615034/understanding-scipy-deconvolve)

