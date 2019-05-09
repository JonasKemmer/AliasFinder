import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.projections import get_projection_class
from matplotlib.transforms import blended_transform_factory
from matplotlib.ticker import MaxNLocator, FixedLocator, NullLocator
from scipy.stats import circmean, circstd


def plot_panel_row(fig, gs, plot_row=0, panel_borders=None,
                   sim_freq=None, gls_obs=None, gls_sim_powers=None,
                   conf_intervals=[[25., 75.], [10., 90.], [2.5, 97.5]],
                   peak_pos=None, obs_phases=None, sim_phases=None,
                   hide_xlabel=True):
    """ Wrapper function to plot one row that contains one simulation.

    Parameters
    ----------
    fig : object
        The figure object from the main script.
    gs : object
        The gridspec from the main script.
    plot_row : int, optional
        Which row from the main figure to plot.
    panel_borders : array, optional
        The borders of the panels showing the different subsets of the GLS.
    sim_freq : float, optional
        The frequency of the simulated signal. Is used to mark it in the plot.
    gls_obs : object, optional
        The GLS of the observed data which is plotted for comparison.
    gls_sim_powers : array, optional
        The array which contains the power information from the simulated GLS.
    conf_intervals : array
        An array which contains the limits of the upper an lower confidence
        intervals which are to be grey shaded around the median of the
        simulated powers. Should be in percentages.
    peak_pos : array
        Array that contains the positions of the detected peaks in the GLS of
        the observed RVs.
    obs_phases : array
        Array containing the phase information about the detected peaks in
        the GLS of the observed RVs.
    sim_phases : array
        Array with the phase information for the peaks in the simulated data.

    """
    axes_row = []
    for panel_idx, panel in enumerate(panel_borders):
        ax = fig.add_subplot(gs[plot_row, panel_idx])
        ax.set_xlim(panel)
        # ax.set_ylim(0, np.max((np.nanmax(np.percentile(gls_sim_powers, 97.5,
        #                                                axis=0)),
        #                        np.nanmax(gls_obs.power)))*1.45)
        if panel_idx == 0:
            ax = plot_info(ax, sim_freq, label=True)
        else:
            ax.set_yticklabels([])
            ax = plot_info(ax, sim_freq)
        ax = plot_lines(ax, gls_obs, gls_sim_powers, intervals=conf_intervals)
        ax = plot_phase_information(ax, peak_pos, obs_phases, sim_phases)
        if np.logical_and(hide_xlabel, plot_row == 2) or not hide_xlabel:
            if len(panel_borders) == 3:
                ax.xaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper',
                                                       min_n_ticks=3))
            else:
                loc = FixedLocator(np.round([panel[0]+0.8*(panel[1]-panel[0])/3,
                                             panel[0]+2.2*(panel[1]-panel[0])/3],
                                            4))
                ax.xaxis.set_major_locator(loc)
        else:
            ax.xaxis.set_major_locator(NullLocator())

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
        axes_row.append(ax)
    return axes_row


def plot_lines(ax, gls_obs, gls_sim_powers,
               intervals=[[25., 75.], [10., 90.], [2.5, 97.5]],
               alph_int=[0.8, 0.5, 0.2]):
    """ Plots the observed data and the simulated data and
        their confidence intervals.

    Parameters
    ----------
    For a description of the parameters see 'plot_panel_row'.

    """
    lower_freq, upper_freq = ax.get_xlim()
    panel_mask = np.logical_and(gls_obs.freq > lower_freq,
                                gls_obs.freq < upper_freq)
    freq = gls_obs.freq[panel_mask]
    sim_powers = np.median(gls_sim_powers, axis=0)[panel_mask]
    ax.plot(freq, sim_powers, color='black', linestyle='-')
    for alpha, percentile in zip(alph_int, intervals):
        lower, upper = np.percentile(gls_sim_powers, percentile,
                                     axis=0)
        ax.fill_between(freq, lower[panel_mask], upper[panel_mask],
                        facecolor='gray', alpha=alpha)
    ax.plot(gls_obs.freq, gls_obs.power, color='red', linestyle='-')
    return ax


def get_theta(phase):
    """ Calculate the theta angle from radians """
    return (phase % 1) * 2 * np.pi


def plot_phase_clock(theta_obs, theta_sim, range_sims, freq,
                     main_axis, width=0.3):
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
    """ Plots the phase information for the detected peaks

    Parameters
    ----------
    For a description of the parameters see 'plot_panel_row'.

    """
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
    """ Plots the marker for the simulated frequency

    Parameters
    ----------
    For a description of the parameters see 'plot_panel_row'.

    """
    # if label:
    #     ax.text(freq+0.2, 0.85, f'P$_\mathrm{{sim.}}$ = {1/freq:.2f} d',
    #             transform=ax.transAxes,
    #             horizontalalignment='left', fontsize='smaller')
    ax.axvline(freq, color='blue', linestyle='--')
    return ax
