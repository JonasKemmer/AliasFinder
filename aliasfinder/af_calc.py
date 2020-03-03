import numpy as np
from scipy.signal import find_peaks

import aliasfinder.af_utils as af_utils


def get_phase_info(gls_obj, power_threshold=None, sim=False,
                   frequency_array=None):
    """ Get information on the phase of all peaks above a given power
        threshold.

    Parameters
    ----------
    gls_obj : object
        A GLS object from Zechmeisters gls.py.
    power_threshold : float, optional
        The power threshold (ZK) above which the peaks should be analyzed.
        Default is 0.5.
    sim : boolean, optional
        If the commited gls object is a simulated one or not.
        Default is False, if it is set to true the frequency array must be
        given too.
    frequency_array : array, optional
        If a simulated periodogram is commited, this array must contain the
        frequencies of the peaks in the GLS of the observed data.

    Returns
    -------
    array
        An array containing the freqencies, powers and phases of all
        analyzed peaks.

    """
    if sim:
        indexes = np.searchsorted(gls_obj.freq, frequency_array)
    else:
        indexes, _ = find_peaks(gls_obj.p, height=power_threshold)
    phases = np.zeros(len(indexes))
    freqs = np.zeros(len(indexes))
    powers = np.zeros(len(indexes))
    for i in range(len(indexes)):
        k = indexes[i]
        pmax = gls_obj.p[k]
        # Best parameters
        fbest = gls_obj.freq[k]
        ph = np.arctan2(gls_obj._a[k], gls_obj._b[k])/(2.*np.pi)
        phases[i] = ph
        freqs[i] = fbest
        powers[i] = pmax
    return [freqs, powers, phases]


def sim_any_sinmode(gls_obj, in_freq, times):
    """ Simulate a sinusoidal signal with a given frequency and amplitude.

    Parameters
    ----------
    gls_obj : object
        A GLS object from Zechmeisters gls.py.
    in_freq : float
        Frequency of the sinus which is calculated.
    times : array
        Time array for which the amplitude of the signal is calculated.

    Returns
    -------
    array
        Array of the amplitudes of the simulated sinusoidal at the input times.

    """
    k = np.where(np.round(gls_obj.freq, 5) == in_freq)[0]
    if len(k) > 1:
        k = k[0]
    amp = np.sqrt(gls_obj._a[k]**2 + gls_obj._b[k]**2)
    ph = np.arctan2(gls_obj._a[k], gls_obj._b[k])/(2.*np.pi)
    T0 = times.min() - ph/in_freq
    offset = gls_obj._off[k] + gls_obj._Y
    return amp * np.sin(2*np.pi*in_freq*(times-T0)) + offset


def sample_gls(times, gls_obs, freq, jitter, rvs_err, fbeg, fend, object_name,
               peaks_data, search_phase_range):
    """The function to calculate a simulated GLS periodogram.

    Parameters
    ----------
    times : array
        The times of observation of the observed RVs.
    gls_obj : object
        A GLS object from Zechmeisters gls.py from the observed RVs.
    freq : float
        The frequency for which the GLS should be simulated.
    jitter : float
        The jitter which is added to the simulated data.
    rvs_err : array
        The errors of the measured RVs. They are adopted as the errors of
        the simulated RVs
    fbeg : float, optional
        Start frequency of the GLS periodogram.
    fend : float, optional
        End frequency of the GLS periodogram.
    object_name : str, optional
        Name of the object from which the RV data originate.
    peaks_data : array
        The array that contains the information about the peaks measured in the
        observed RVs. The frequencies are used to get information about these
        peaks in the simulated GLS
    search_phase_range : float
        The range around the peaks in the observed data for which peaks in
        the simulated GLS are searched for.

    Returns
    -------
    arrays
        Returns the powers at the corresponding frequencies of the simulated
        GLS as well as the phases of the peaks in the simulated GLS.

    """
    np.random.seed()
    rvs_sim = np.ones(np.shape(times)) \
              + sim_any_sinmode(gls_obs, freq, times) \
              + np.random.normal(0, np.sqrt(jitter**2), times.size)
    ls_sim = af_utils.get_gls(times, rvs_sim, rvs_err, fbeg, fend, object_name,
                              freq_array=gls_obs.freq)
    dummy_freq_array = np.zeros(np.size(peaks_data[0]))
    # search for phases of max power using a certian frequency
    # range and the prior of data peaks
    for j in range(0, np.size(peaks_data[0])):
        index_frequencies = np.where(np.logical_and(ls_sim.freq >= peaks_data[0][j]
                                                    - search_phase_range,
                                                    ls_sim.freq <= peaks_data[0][j]
                                                    + search_phase_range))
        index_maxfreqs = max(np.arange(len(ls_sim.p[index_frequencies])),
                             key=ls_sim.p[index_frequencies].__getitem__)
        index_maxfreq = np.argwhere(ls_sim.freq ==
                                    ls_sim.freq[index_frequencies][index_maxfreqs])[0]

        dummy_freq_array[j] = ls_sim.freq[index_maxfreq]
    peaks_sim = get_phase_info(ls_sim,
                                       frequency_array=dummy_freq_array,
                                       sim=True)
    return ls_sim.power, ls_sim.freq, (peaks_sim[2] % 1) * 2. * np.pi

def get_metric(gls_obs=None, gls_sim_powers=None):
    """ The function to calculate the likelihood.

    Parameters
    ----------
    gls_obs : object, optional
        The GLS of the observed data.
    gls_sim_powers : array, optional
        The array which contains the power information from the simulated GLS.
    conf_intervals : array, optional
        An array which contains the limits of the upper an lower confidence
        intervals which are to be grey shaded around the median of the
        simulated powers. Should be in percentages.
    peak_pos : array, optional
        Array that contains the positions of the detected peaks in the GLS of
        the observed RVs.
    obs_phases : array, optional
        Array containing the phase information about the detected peaks in
        the GLS of the observed RVs.
    sim_phases : array, optional
        Array with the phase information for the peaks in the simulated data.

    """

    freq = gls_obs.freq[:]
    data_powers = gls_obs.power[:]
    sim_powers = np.median(gls_sim_powers, axis=0)
    sigma_powers = np.std(gls_sim_powers,axis=0)

    #Li=np.exp(-(data_powers-sim_powers)**2/sigma_powers**2,dtype=np.float128)
    #chisq=(data_powers-sim_powers)**2/sigma_powers**2
    squarediff=(data_powers-sim_powers)**2
    L=np.sum(squarediff)


    return L
