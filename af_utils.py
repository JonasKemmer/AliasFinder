import numpy as np
import GLS
from detect_peaks import detect_peaks
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes #for the plotting of polar axes in plots
from matplotlib.projections import get_projection_class
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.transforms import blended_transform_factory
from matplotlib.ticker import MaxNLocator, IndexLocator, MultipleLocator
from scipy.stats import circmean, circstd


def find_nearest(array, value):
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
    """ Reads the RV data into different arrays"""
    alldata = []
    for k, f in enumerate(rv_files):
        alldata.append(np.genfromtxt("/%s"%(f), skip_header=0, unpack=True,
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
    """ Wrapper to call GLS.py from Mathias Zechmeister """
    if freq_array != []:
        ls = GLS.Gls((t, rv, err), verbose=verbose, fbeg=fbeg, fend=fend,
                     name=object_name, norm=norm, ls=ls, ofac=20,
                     freq=freq_array)
    else:
        ls = GLS.Gls((t, rv, err), verbose=verbose, fbeg=fbeg, fend=fend,
                     name=object_name, norm=norm, ls=ls, ofac=20)

    return ls


def get_spectral_window_function(freqs, times):
    """ Function to calculate the spectral window function of an series
        of observations given an array fo test frequencies."""
    w_v = np.zeros(freqs.shape)
    for idx, freq in enumerate(freqs):
        w_v[idx] = np.sum(np.exp(-2j*np.pi*freq*times))/len(times)
    return np.abs(w_v)


def get_phase_info(gls_obj, frequency_array=None, power_threshold=None,
                   sim=True):
    """ get infomartion on the phase of all peaks above a given power
        threshold"""
    if sim:
        perio_peaks = peaksPeriodogram(gls_obj, power_threshold)
    else:
        perio_peaks = peaksSimulated(gls_obj, frequency_array)
    return perio_peaks


def peaksPeriodogram(GLS, power_threshold):
    """
    Analyze the highest periodogram peaks with a FAP over a certain threshold.
    """
    indexes = detect_peaks(GLS.p, mph=power_threshold)
    phases = np.zeros(len(indexes))
    freqs = np.zeros(len(indexes))
    powers = np.zeros(len(indexes))
    for i in range(len(indexes)):
        k = indexes[i]
        pmax = GLS.p[k]

        # Best parameters
        fbest = GLS.freq[k]
        ph = np.arctan2(GLS._a[k], GLS._b[k])/(2.*np.pi)
        phases[i] = ph
        freqs[i] = fbest
        powers[i] = pmax
        if 1 < k < GLS.nf-2:
            # Shift the parabola origin to power peak
            xh = (GLS.freq[k-1:k+2] - GLS.freq[k])**2
            yh = GLS.p[k-1:k+2] - pmax
            # Calculate the curvature (final equation from least square)
            aa = np.dot(yh, xh) / np.dot(xh, xh)
        else:
            print('WARNING: Highest peak is at the edge of the frequency'
                  ' range.\nNo output of frequency error.\nIncrease '
                  'frequency.')
    return [freqs, powers, phases]


def peaksSimulated(GLS, frequency_array):
    """
    Analyze the highest periodogram peaks with a FAP over a certain threshold.
    """
    indexes = np.searchsorted(GLS.freq, frequency_array)
    phases = np.zeros(len(indexes))
    freqs = np.zeros(len(indexes))
    powers = np.zeros(len(indexes))
    for i in range(len(indexes)):
        k = indexes[i]
        # Power
        pmax = GLS.p[k]
        # Best parameters
        fbest = GLS.freq[k]
        ph = np.arctan2(GLS._a[k], GLS._b[k])/(2.*np.pi)
        phases[i] = ph
        freqs[i] = fbest
        powers[i] = pmax

        # Get the curvature in the power peak by fitting a parabola y=aa*x^2
        if 1 < k < GLS.nf-2:
            # Shift the parabola origin to power peak
            xh = (GLS.freq[k-1:k+2] - GLS.freq[k])**2
            yh = GLS.p[k-1:k+2] - pmax
            # Calculate the curvature (final equation from least square)
            aa = np.dot(yh, xh) / np.dot(xh, xh)

        else:
            print('WARNING: Highest peak is at the edge of the frequency'
                  ' range.\nNo output of frequency error.\nIncrease '
                  'frequency.')
    return [freqs, powers, phases]


def sim_any_sinmode(gls, in_freq, times):
    """ Simulate a np.sinusoidal signal with a given frequency and
        amplitude """
    k = np.where(np.round(gls.freq, 5) == in_freq)[0]
    amp = np.sqrt(gls._a[k]**2 + gls._b[k]**2)
    ph = np.arctan2(gls._a[k], gls._b[k])/(2.*np.pi)
    T0 = times.min() - ph/in_freq
    offset = gls._off[k] + gls._Y
    return amp * np.sin(2*np.pi*in_freq*(times-T0)) + offset


def plot_panel_row(fig, gs, plot_row=0, panel_borders=None,
                   sim_freq=None, gls_obs=None, gls_sim_powers=None,
                   conf_intervals=[[25., 75.], [10., 90.], [2.5, 97.5]],
                   peak_pos=None, obs_phases=None, sim_phases=None):
    """ Wrapper function to plot the one row that contains one simulation """
    for panel_idx, panel in enumerate(panel_borders):
        ax = fig.add_subplot(gs[plot_row, panel_idx])
        ax.set_xlim(panel)
        ax.set_ylim(0, np.max((np.nanmax(gls_sim_powers),
                               np.nanmax(gls_obs.p)))*1.4)
        if panel_idx == 0:
            ax = plot_info(ax, sim_freq, label=True)
        else:
            ax.set_yticklabels([])
            ax = plot_info(ax, sim_freq)
        ax = plot_lines(ax, gls_obs, gls_sim_powers, intervals=conf_intervals)
        ax = plot_phase_information(ax, peak_pos, obs_phases, sim_phases)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper',
                                               min_n_ticks=3))
        d = .015
        if panel_idx == 0:
            ax.tick_params(labelright='off')
            kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
            ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
            ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
            ax.spines['right'].set_visible(False)
        elif panel_idx < len(panel_borders)-1:
            ax.set_yticks([])
            kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
            ax.plot((-d, +d), (-d, +d), **kwargs)
            ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)
            ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
            ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
        else:
            ax.set_yticks([])
            ax.tick_params(labelright='off')
            ax.tick_params(labelleft='off')
            kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
            ax.plot((-d, +d), (-d, +d), **kwargs)
            ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)
            ax.spines['left'].set_visible(False)



def plot_lines(ax, gls_obs, gls_sim_powers,
               intervals=[[25., 75.], [10., 90.], [2.5, 97.5]]):
    """ Plots the observed data and the simulated data and
        their confidence intervals """
    lower_freq, upper_freq = ax.get_xlim()
    panel_mask = np.logical_and(gls_obs.freq > lower_freq,
                                gls_obs.freq < upper_freq)
    freq = gls_obs.freq[panel_mask]
    sim_powers = np.median(gls_sim_powers, axis=0)[panel_mask]
    ax.plot(freq, sim_powers, color='black', linestyle='-')
    for percentile in intervals:
        lower, upper = np.percentile(gls_sim_powers, percentile,
                                     axis=0)
        ax.fill_between(freq, lower[panel_mask], upper[panel_mask],
                        facecolor='gray', alpha=0.5)
    ax.plot(gls_obs.freq, gls_obs.power, color='red', linestyle='-')
    return ax


def get_theta(phase):
    """ Calculate the theta angle from radians """
    return (phase % 1) * 2 * np.pi


def plot_phase_clock(theta_obs, theta_sim, range_sims, freq,
                     main_axis, width=0.25):
    """ Plot function for an inset phase-clock """
    trans = blended_transform_factory(main_axis.transData,
                                      main_axis.transAxes)
    ax_sub = inset_axes(main_axis, width=width, height=width, loc=10,
                        bbox_to_anchor=(freq, 0.85),
                        bbox_transform=trans,
                        borderpad=0.0,
                        axes_class=get_projection_class("polar"))
    radii = [90]
    ax_sub.vlines(theta_obs, 0, radii, color='red')
    ax_sub.vlines(theta_sim, 0, radii, color='black')
    ax_sub.fill_between(np.linspace(range_sims[0], range_sims[1], 100), 0,
                        radii, color='gray', alpha=0.6)
    ax_sub.set_theta_direction(-1)
    ax_sub.set_theta_offset(np.pi/2.0)
    ax_sub.set_thetagrids([])
    ax_sub.set_yticks([])
    return main_axis


def plot_phase_information(ax, peak_pos, obs_phases, sim_phases):
    """ Plots the phase information for the detected peaks"""
    sim_mean_phases = circmean(sim_phases, axis=0)
    sim_phase_stds = circstd(sim_phases, axis=0)
    sim_phase_ranges = np.array((sim_mean_phases-sim_phase_stds,
                                 sim_mean_phases+sim_phase_stds)).T
    lower_freq, upper_freq = ax.get_xlim()
    for freq, obs_phase, sim_phase, sim_phase_range in zip(peak_pos,
                                                           obs_phases,
                                                           sim_mean_phases,
                                                           sim_phase_ranges):
        if lower_freq < freq < upper_freq:
            plot_phase_clock(get_theta(obs_phase), sim_phase,
                             sim_phase_range,
                             freq, ax)
    return ax


def plot_info(ax, freq, label=False):
    """ Plots the marker for the simulated frequency"""
#    if label:
#        ax.text(0.05, 1.05, f'sim. P = {1/freq:.4f} d', transform=ax.transAxes,
#                horizontalalignment='left', fontsize='smaller')
    ax.axvline(freq, color='grey', linestyle='--')
    return ax
