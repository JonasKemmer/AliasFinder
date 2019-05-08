import os
import sys
import warnings
from multiprocessing import Pool

import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.widgets import Cursor
from tqdm import tqdm

import af_calc
import af_plots
import af_utils


__author__ = "Jonas Kemmer, Stephan Stock @ ZAH, Landessternwarte Heidelberg"
__version__ = "1.0."
__license__ = "MIT"


warnings.filterwarnings("ignore")

if __name__ == '__main__':
    print(f' \n'
          f'Welcome to the Alias Finder V{__version__}\n')
    # Input the parameter set in the parameter file - the list of possible
    # arguments is defined in the default_params-dictionairy
    with open(str(sys.argv[1])) as stream:
        try:
            input_params = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    default_params = {'object_name': '',
                      'savepath': './',
                      'rv_files': [],
                      'test_period': 50,
                      'sampling_freq': None,
                      'mc_samples': 100,
                      'num_processes': 1,
                      'alias_order': 1,
                      'panel_width': 1,
                      'hide_xlabel': True,
                      'fbeg': 0.0001,
                      'fend': 2.5,
                      'power_threshold': 0.05,
                      'search_phase_range': 0.00005,
                      'offsets': [0],
                      'substract': False,
                      'use_rms_as_jitter': True,
                      'jitter': 0}

    for key in default_params:
        try:
            locals().update({key: input_params[key]})
        except KeyError:
            locals().update({key: default_params[key]})
            print(f'Parameter {key} is missing in the input file.\n'
                  f'Value set to default: {key} = {default_params[key]}\n')

    # Read in the observed RV data
    times, rvs, rvs_err = af_utils.read_rvs(rv_files, offsets)
    if sampling_freq is None or sampling_freq == 'None':
        # First analyze the spectral window function of the data and select a
        # sampling frequency who's aliases should be analyzed
        window_freqs = np.linspace(fbeg, fend, len(times)*30)
        window_powers = af_utils.get_spectral_window_function(window_freqs,
                                                              times)

        print('Please select the sampling frequency you want to test '
              'from the plot')
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plt.title('Spectral window function of the data\n'
                  'Please select the sampling frequency you want to test '
                  '(Press any button then click)')
        plt.xlabel('Frequency f [1/d]')
        plt.ylabel(r'Power W($\nu$)')
        ax.plot(window_freqs, window_powers, color='black')
        cursor = Cursor(ax, useblit=True, color='k', linewidth=1)
        zoom_ok = False
        print('\nZoom or pan to view, \npress any button when ready'
              ' to select the sampling frequency - then click:\n')
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress()
        sampling_freq = plt.ginput(1)
        sampling_freq = np.round(sampling_freq[0][0], 5)
        plt.axis('auto')
        plt.autoscale()
        plt.xlim(0, 1.2)
        plt.title('')
        fig.set_size_inches(3.5, 3, forward=True)
        axins = ax.inset_axes([0.2, 0.53, 0.4, 0.4])
        axins.plot(window_freqs, window_powers, color='black')
        axins.set_xlim(0.99, 1.01)
        axins.xaxis.label.set_fontsize(4)
        axins.yaxis.label.set_fontsize(4)
        plt.tight_layout()
        save_var = os.path.join(savepath, f'{object_name}_spectral_window_'
                                          'function.pdf')
        plt.savefig(save_var,
                    bbox_inches='tight', dpi=400)
        plt.close()
        print(f'Your selected sampling frequency is {sampling_freq}')
    else:
        print(f'Sampling frequency (f = {sampling_freq}) given in the input '
              'file.\nSkipping selection from plot of the window function.\n'
              'If you want to select the sampling frequeny from the window '
              'function\nplease set the it to "None" in the input file.\n')

    # Second, create an GLS from the observed RVs
    gls_obs = af_utils.get_gls(times, rvs, rvs_err, fbeg, fend,
                               object_name)

    # If a Kepler signal is to be substracted:
    if substract == 1:
        print('You chose to substract the most significant Period!')
        rv_array = rvs - gls_obs.sinmod(times)
        gls_obs = af_utils.get_gls(times, rv_array, rvs_err, fbeg, fend,
                                   object_name, verbose=False)
        rvs = rv_array
        print('Above Signal now significant!')

    if substract == 2:
        print('You chose to substract the 2 most significant period!')
        rv_array = rvs - gls_obs.sinmod(times)
        gls_obs = af_utils.get_gls(times, rv_array, rvs_err, fbeg, fend,
                                   object_name, verbose=False)
        rv_array = rv_array-gls_obs.sinmod(times)
        gls_obs = af_utils.get_gls(times, rv_array, rvs_err, fbeg, fend,
                                   object_name, verbose=False)
        rvs = rv_array
        print('Above Signal now most significant!')

    # Find peaks in the GLS of the observed RVs
    peaks_data = af_calc.get_phase_info(gls_obs,
                                        power_threshold=power_threshold)

    # Select the Periods which are to be simulated
    sim_freqs = np.zeros(3)
    for idx, m in enumerate([0, -1, 1]):
        while True:
            wanted_freq = np.round(af_utils.calc_alias_freq(test_period,
                                                            sampling_freq, m),
                                   5)
            closest_freq = np.round(af_utils.find_nearest_frequency(peaks_data[0],
                                                                    wanted_freq),
                                    5)
            check = input(f'Please check {idx+1}. frequency you want to simulate\n'
                          f'The predicted freq is:            {wanted_freq}\n'
                          f'The closest found freq is:        {closest_freq}\n'
                          f'Do they fit? (yes/no)\n'
                          f'(If you want to manually choose another \n'
                          f'frequency please also type \'no\')  \n')
            if check.lower() == 'yes' or check.lower() == 'y':
                sim_freqs[idx] = closest_freq
                break
            elif check.lower() == 'no' or check.lower() == 'n':
                print('The alias was not found in the data. Try to lower \n'
                      'the power threshold to find more peaks.')
                by_hand = input('Do you want to choose manually? (yes/no)')
                if by_hand.lower() == 'yes' or by_hand.lower() == 'y':
                    for dim, key in zip(peaks_data, ['Frequencies', 'Powers',
                                                     'Phases']):
                        print(f'{key}\n {np.round(dim, 5)}\n')
                    sim_freqs[idx] = input('The first array shows the frequencies'
                                           ', the second the power, the third the '
                                           'phases. Please chose the frequency '
                                           'that is closest to the one you want '
                                           'to simulate (as float, rounded to the '
                                           'first 5 digits, in 1/days): ')
                    if not isinstance(sim_freqs[idx], float):
                        print("You need to specify it as float! Try again.")
                    if np.isin(sim_freqs[idx], np.round(peaks_data[0], 5)):
                        print(f'You chose f = {sim_freqs[idx]}\n\n')
                        break
                    else:
                        print('Your chosen frequency is not in the peak array.'
                              ' Try again.')
                else:
                    raise Exception('No frequency to simulate!')

    # Definition of the borders of the plot panels
    f1s = 1/(test_period+panel_width)
    f1e = 1/(test_period-panel_width)
    f2s = np.abs(1/(test_period-panel_width)-sampling_freq)
    f2e = np.abs(1/(test_period+panel_width)-sampling_freq)
    f3s = np.abs(1/(test_period+panel_width)+sampling_freq)
    f3e = np.abs(1/(test_period-panel_width)+sampling_freq)
    if alias_order == 1:
        panels = np.array([[f1s, f1e], [f2s, f2e], [f3s, f3e]])
    elif alias_order == 2:
        panels = np.array([[f1s, f1e], [f2s, f2e], [f3s, f3e], [f2s+1, f2e+1],
                           [f3s+1, f3e+1]])
    else:
        raise Exception('Orders higher than m=2 are not implemented.')

    # Initializing the plot and starting the simulations
    fig = plt.figure(figsize=(7, 4))
    if alias_order == 1:
        gs = gridspec.GridSpec(3, 3)
    else:
        gs = gridspec.GridSpec(3, 5)

    if hide_xlabel:
        gs.update(wspace=0.08, hspace=0)
    else:
        gs.update(wspace=0.08, hspace=0.25)

    matplotlib.rc('xtick', labelsize=9)
    matplotlib.rc('ytick', labelsize=9)
    fig.text(0.5, 0.02, 'Frequency f [1/d]', ha='center')
    fig.text(0.05, 0.5, gls_obs.label["ylabel"], va='center',
             rotation='vertical')

    if use_rms_as_jitter:
        jitter = gls_obs.rms
    else:
        jitter = jitter

    for sim_freq_idx, freq in enumerate(sim_freqs):
        print(f'\n Simulating Freq. # {sim_freq_idx+1:.0f}\n')
        gls_sim_powers = np.zeros((mc_samples, len(gls_obs.power)))
        freqs_array = np.zeros((mc_samples, len(gls_obs.freq)))
        sim_phases = np.zeros((mc_samples, len(peaks_data[2])))
        if num_processes > 1:
            pbar = tqdm(total=mc_samples)
            def update(*a):
                pbar.update()
            with Pool(num_processes) as p:
                res_temp = [p.apply_async(af_calc.sample_gls,
                                          args=(times, gls_obs, freq,  jitter,
                                                rvs_err, fbeg,  fend,
                                                object_name, peaks_data,
                                                search_phase_range),
                                          callback=update)
                            for i in range(pbar.total)]
                for idx, result in enumerate(res_temp):
                    gls_sim_powers[idx], \
                     freqs_array[idx], \
                     sim_phases[idx] = result.get()
        else:
            for i in tqdm(range(0, mc_samples)):
                gls_sim_powers[i], \
                 freqs_array[i], \
                 sim_phases[i] = af_calc.sample_gls(times, gls_obs, freq,
                                                    jitter, rvs_err, fbeg,
                                                    fend, object_name,
                                                    peaks_data,
                                                    search_phase_range)
        af_plots.plot_panel_row(fig, gs, plot_row=sim_freq_idx,
                                panel_borders=panels, sim_freq=freq,
                                gls_obs=gls_obs, gls_sim_powers=gls_sim_powers,
                                peak_pos=peaks_data[0],
                                obs_phases=peaks_data[2],
                                sim_phases=sim_phases,
                                hide_xlabel=hide_xlabel)

    plt.tight_layout()
    save_var = os.path.join(savepath, f'{object_name}_{test_period}d'
                            '_alias_test.pdf')
    plt.savefig(save_var,
                bbox_inches='tight', dpi=400)
    plt.show()
    print(f'\n \n Plot saved to {save_var}')
