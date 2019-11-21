---
title: 'AliasFinder: A Python script to search for the true planetary frequency within radial velocity data'
tags:
  - Python
  - astronomy
  - exoplanets
  - radial velocity
  - signal search
authors:
  - name: Stephan Stock
    orcid: 0000-0002-1166-9338
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Jonas Kemmer
    orcid: 0000-0003-3929-1442
    affiliation: "1, 2"
affiliations:
 - name: Landessternwarte, Zentrum für Astronomie der Universität Heidelberg, Königstuhl 12, 69117 Heidelberg, Germany
   index: 1
 - name: Fellow of the International Max Planck Research School for Astronomy and Cosmic Physics at the University of Heidelberg (IMPRS-HD)
   index: 2
date: 23 July 2019
bibliography: paper_alias.bib
---

# Summary

Planets around other stars, so called exoplanets, can be detected with the radial
velocity (RV) method which uses high-resolution spectroscopy on the host star. The method transforms the
measured Doppler-shift in a stellar spectrum to a radial velocity.
This Doppler-shift is caused by the stellar reflex motion (wobbling),
in particular the line of sight movement of the star, due to its gravitational interaction with the orbiting planet.
By now the RV method is the second most successful method and the first one to find exoplanets.
The detection of exoplanets by this method relies
on finding periodic signals in noisy time-series data (RVs over time) using a periodogram analysis.
However, gaps in the observations (e.g. stars can only be observed during night from the ground)
can result into strong aliases, an effect where two signals are almost indistinguishable when sampled. Aliases
normally occure at $f_{\text{Alias}}=f_{\text{True}}\pm m \cdot f_{\text{Sampling}}$ where $f_{\text{True}}$
is the true signal, $f_{\text{Sampling}}$ the sampling frequency and $m$ an integer value.
This makes it difficult to asses the true period of the planet and has already led to incorrectly published
orbital periods (see Dawson and Fabrcky 2010; Stock et al. in prep).

``AliasFinder`` is a Python script which can help to identify the true signal within a time series data. 
The script is based on a method described in @Dawson2010, however ``AliasFinder`` uses some improvements described
in detail within Stock et al.(in prep.). **To our knowledge, there are no other public packages to perform alias testing based on this methodology.** The basic principle is to simulate an ensemble of noisy time-series based on the
original data sampling and noise level and inject one of the signals observed within the original data into this time-series.
 ``AliasFinder`` will than plot the median GLS periodogram [@Zechmeister2009] from the simulations of the injected
signal overplotted with the original data periodogram. From the comparison of the signal properties (signal power,
frequeny or phase) at the injected frequency and its aliases between these two periodograms one can
than asses the true period of the signal. We show a straightforward example of such a plot in Fig. 1. where we test two possible true periods of a
planet at either 2.64 d or 1.60 d which was analzed by [@Trifonov2018]. From the comparison of the simulated time series (in black) to the observed periodogram (in red) it is obvious that only the period of 2.64 d (f~0.38 1/d) is able to reproduce the data correctly.

``AliasFinder`` is specifically written for Python 3 and should run with a standard Anaconda distribution
with a few additional packages. The script is executed by passing a yaml parameter file with important information
about the star, planetary system and the measured RV signals.  

``AliasFinder`` is developed to be used by astronomers that search for exoplanets in RV data. Nevertheless
the software can in principle be used by anyone who needs to distinguish a possible signal from its aliases
within any noisy time-series if the correct input format is respected. The script is already
used within Stock et al. (in prep.) and will be used in further studies in the future.
``AliasFinder`` is designed to be fast (multi-core support), user-friendly and produce high quality publication plots as in Fig. 1.
This will encourage the usage of this test within future scientific publications based on time-series data,
in particular RV data, and result in less incorrectly published orbital periods for exoplanets.

![Example of a test for aliases. The top panel shows the resulting median periodogram (in black) derived from simulations of an ensemble of GLS periodograms where a signal of 2.64 d (marked with the vertical blue dashed line) was injected which was also observed within the original data periodogram (in red).
The two lower panels show the simulated median periodogram where the first order aliases of 2.64 d which are at 1.6 d and 0.72 d are injected (also marked with the blue dashed line). The periodograms in each row are automatically plotted to regions that allow to distinguish the true period from its first oder aliases. The interquartile range and the range of 90\% and 99\% of the simulated periodograms are shown as the shaded grey area.
The angular mean of the phase of each peak and its standard deviation
are shown in the clock diagram (red: simulated; grey: data).](example2.png)

# Acknowledgements]

We acknowledge contributions from Paul Heeren, and support from Sabine Reffert,
and Trifon Trifonov during the genesis of this project. **We made us of NumPy [@oliphant2006guide]. This research made use of Astropy [@Astropy_comm] a community-developed core Python package for Astronomy [@Astropy2018]. We made use of tqdm [@tqdm], SciPy @[Scipy2019] and matplotlib [@Hunter2007].**

# References
