import numpy as np

import aliasfinder.af_gls as af_gls


def find_nearest_frequency(array, value):
    """ Helper function to find the entry in a array closest to the given
        value"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def calc_alias_freq(freq, freq_s, m):
    """ Calculates the alias given the test frequency,
        the sampling frequency and the order m."""
    return np.abs(1/(freq)+m*freq_s)


def read_rvs(rv_files, offsets):
    """ Reads the RV data into different arrays.

    Parameters
    ----------
    rv_files : list
        A list of strings containing the paths to the RV files.
        First column must be the times, second the RV and third the RV errors
    offsets : list
        If multiple datasets are given this list should contain the offsets
        between the datasets.

    Returns
    -------
    np.arrays
        Returns three numpy array containing the times, rvs and rv-errors.

    """

    alldata = []
    for k, f in enumerate(rv_files):
        alldata.append(np.genfromtxt(f, skip_header=0, unpack=True,
                       skip_footer=0))
    times = np.hstack(alldata[j][0] for j in range(len(alldata)))
    rvs = np.hstack(alldata[j][1] - offsets[j] for j in range(len(alldata)))
    errors = np.hstack(alldata[j][2] for j in range(len(alldata)))
    indexes = np.argsort(times)
    times = times[indexes]
    rvs = rvs[indexes]
    errors = errors[indexes]
    return np.array(times), np.array(rvs), np.array(errors)


def get_gls(t, rv, err, fbeg=0.000001, fend=1.06, object_name=None,
            norm='ZK', ls=False, freq_array=[], verbose=False):
    """ Wrapper to call GLS.py from Mathias Zechmeister.
        https://github.com/mzechmeister/GLS/blob/master/python/gls.py

    Parameters
    ----------
    t : array
        Array containing the times of the observations.
    rv : array
        Array containing the Rvs.
    err : array
        Array with the errors of the Rvs.
    fbeg : float, optional
        Start frequency of the GLS periodogram.
    fend : float, optional
        End frequency of the GLS periodogram.
    object_name : str, optional
        Name of the object from which the RV data originate.
    norm : string, optional
        The normalization; either of "ZK", "Scargle", "HorneBaliunas",
        "Cumming", "wrms", "chisq". The default is unity ("ZK").
    ls : boolean, optional
        If True, the conventional Lomb-Scargle periodogram will be computed
        (default is False).
    freq_array : array, optional
        Contains the frequencies at which to calculate the periodogram.
        If given, fast and verbose option are not available.
        If not given, a frequency array will be automatically generated.
    verbose : boolean, optional
        Set True to obtain some statistical output (default is False).

    Returns
    -------
    object
        A GLS object from Matthias Zechhmeisters Script.

    """
    if freq_array != []:
        ls = af_gls.Gls((t, rv, err), verbose=verbose, fbeg=fbeg, fend=fend,
                        name=object_name, norm=norm, ls=ls, ofac=20,
                        freq=freq_array)
    else:
        ls = af_gls.Gls((t, rv, err), verbose=verbose, fbeg=fbeg, fend=fend,
                        name=object_name, norm=norm, ls=ls, ofac=20)

    return ls


def get_spectral_window_function(freqs, times):
    """ Function to calculate the spectral window function of an series
        of observations given an array fo test frequencies.
        (Roberts et al. 1987, Dawson & Fabrycky 2018)"""
    w_v = np.zeros(freqs.shape)
    for idx, freq in enumerate(freqs):
        w_v[idx] = np.sum(np.exp(-2j*np.pi*freq*times))/len(times)
    return np.abs(w_v)
