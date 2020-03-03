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
from matplotlib.ticker import MaxNLocator
from tqdm import tqdm
from astropy.table import Table


import aliasfinder.af_calc as af_calc
import aliasfinder.af_plots as af_plots
import aliasfinder.af_utils as af_utils


__author__ = "Jonas Kemmer, Stephan Stock @ ZAH, Landessternwarte Heidelberg"
__version__ = "1.0"
__license__ = "MIT"


warnings.filterwarnings("ignore")

def main():
    print(f' \n'
          f'Welcome to the AliasFinder V{__version__}\n')
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
                      'plot_additional_period_axis': True,
                      'fbeg': 0.0001,
                      'fend': 2.5,
                      'power_threshold': 0.05,
                      'search_phase_range': 0.00005,
                      'offsets': [0],
                      'substract': False,
                      'use_rms_as_jitter': True,
                      'jitter': 0,
                      'calc_metric': False}
    params = {}
    for key in default_params:
        try:
            params[key] = input_params[key]
        except KeyError:
            params[key] = default_params[key]
            print(f'Parameter {key} is missing in the input file.\n'
                  f'Value set to default: {key} = {default_params[key]}\n')
    # Read in the observed RV data
    times, rvs, rvs_err = af_utils.read_rvs(params['rv_files'],
                                            params['offsets'])
    if params['sampling_freq'] is None or params['sampling_freq'] == 'None':
        # First analyze the spectral window function of the data and select a
        # sampling frequency who's aliases should be analyzed
        window_freqs = np.linspace(params['fbeg'], params['fend'],
                                   len(times)*25)
        window_powers = af_utils.get_spectral_window_function(window_freqs,
                                                              times)

        print('Please select the sampling frequency you want to test '
              'from the plot')
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plt.title('Spectral window function of the data\n'
                  'Please select the sampling frequency you want to test '
                  '(Double click left mouse button)')
        plt.xlabel('Frequency f [1/d]')
        plt.ylabel(r'Power W($\nu$)')
        main_line = ax.plot(window_freqs, window_powers, color='black')
        cursor = Cursor(ax, useblit=True, color='k', linewidth=1)

        def onclick(event):
            if event.button == 1 and event.dblclick:
                params['sampling_freq'] = np.round(event.xdata, 5)
                plt.axis('auto')
                plt.autoscale()
                plt.xlim(0, 1.2)
                plt.title('')
                fig.set_size_inches(3.5, 3, forward=True)
                axins = ax.inset_axes([0.2, 0.53, 0.4, 0.4])
                axins.plot(window_freqs, window_powers, color='black',
                           linewidth=0.6)
                axins.set_xlim(0.99, 1.01)
                plt.setp(main_line, linewidth=0.6)
                ax.xaxis.label.set_fontsize(8)
                ax.yaxis.label.set_fontsize(8)
                axins.tick_params(labelsize=8)
                ax.tick_params(labelsize=8)
                plt.tight_layout()
                save_string = os.path.join(params['savepath'],
                                           params['object_name']+'_spectral_window_'
                                           'function.pdf')
                plt.savefig(save_string,
                            bbox_inches='tight', dpi=400)
                plt.close()

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()

        if params['sampling_freq'] is None or params['sampling_freq'] == 'None':
            raise Exception('Plot closed before sampling frequency was '
                            'choosen')
        print('Your selected sampling frequency is '
              + str(params['sampling_freq']))
    else:
        print('Sampling frequency (f = '+str(params['sampling_freq'])
              + ') given in the input '
              'file.\nSkipping selection from plot of the window function.\n'
              'If you want to select the sampling frequeny from the window '
              'function\nplease set the it to "None" in the input file.\n')

    # Second, create an GLS from the observed RVs
    gls_obs = af_utils.get_gls(times, rvs, rvs_err, params['fbeg'],
                               params['fend'], params['object_name'])

    # If a Kepler signal is to be params['substract']ed:
    if params['substract'] == 1:
        print('You chose to substract the most significant Period!')
        rv_array = rvs - gls_obs.sinmod(times)
        gls_obs = af_utils.get_gls(times, rv_array, rvs_err, params['fbeg'],
                                   params['fend'], params['object_name'],
                                   verbose=False)
        rvs = rv_array
        print('Above Signal now significant!')

    if params['substract'] == 2:
        print('You chose to substract the 2 most significant period!')
        rv_array = rvs - gls_obs.sinmod(times)
        gls_obs = af_utils.get_gls(times, rv_array, rvs_err, params['fbeg'],
                                   params['fend'], params['object_name'],
                                   verbose=False)
        rv_array = rv_array-gls_obs.sinmod(times)
        gls_obs = af_utils.get_gls(times, rv_array, rvs_err, params['fbeg'],
                                   params['fend'], params['object_name'],
                                   verbose=False)
        rvs = rv_array
        print('Above Signal now most significant!')

    # Find peaks in the GLS of the observed RVs
    peaks_data = af_calc.get_phase_info(gls_obs,
                                        power_threshold=params['power_threshold'])

    # Select the Periods which are to be simulated
    sim_freqs = np.zeros(3)
    for idx, m in enumerate([0, -1, 1]):
        while True:
            wanted_freq = np.round(af_utils.calc_alias_freq(params['test_period'],
                                                            params['sampling_freq'], m),
                                   5)
            closest_freq = np.round(af_utils.find_nearest_frequency(peaks_data[0],
                                                                    wanted_freq),
                                    5)
            check = input(f'Please check {idx+1}. frequency you want to simulate\n'
                          f'The predicted freq is:            {wanted_freq} '
                          f'(= {(1/wanted_freq):.4f} d)\n'
                          f'The closest found freq is:        {closest_freq} '
                          f'(= {(1/closest_freq):.4} d)\n'
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
    f1s = 1/(params['test_period']) - params['panel_width']
    f1e = 1/(params['test_period']) + params['panel_width']
    f2s = np.abs(1/(params['test_period'])-params['sampling_freq']) - params['panel_width']
    f2e = np.abs(1/(params['test_period'])-params['sampling_freq']) + params['panel_width']
    f3s = np.abs(1/(params['test_period'])+params['sampling_freq']) - params['panel_width']
    f3e = np.abs(1/(params['test_period'])+params['sampling_freq']) + params['panel_width']
    # f1s = 1/(params['test_period']+params['panel_width'])
    # f1e = 1/(params['test_period']-params['panel_width'])
    # f2s = np.abs(1/(params['test_period']-params['panel_width'])-params['sampling_freq'])
    # f2e = np.abs(1/(params['test_period']+params['panel_width'])-params['sampling_freq'])
    # f3s = np.abs(1/(params['test_period']+params['panel_width'])+params['sampling_freq'])
    # f3e = np.abs(1/(params['test_period']-params['panel_width'])+params['sampling_freq'])
    if params['alias_order'] == 1:
        panels = np.array([[f1s, f1e], [f2s, f2e], [f3s, f3e]])
    elif params['alias_order'] == 2:
        panels = np.array([[f1s, f1e], [f2s, f2e], [f3s, f3e], [f2s+1, f2e+1],
                           [f3s+1, f3e+1]])
    else:
        raise Exception('Orders higher than m=2 are not implemented.')

    # Initializing the plot and starting the simulations
    fig = plt.figure(figsize=(7.3, 3.1))
    plt.rcParams.update({'lines.linewidth':1})
    if params['alias_order'] == 1:
        gs = gridspec.GridSpec(3, 3)
    else:
        gs = gridspec.GridSpec(3, 5)

    if params['hide_xlabel']:
        gs.update(wspace=0.08, hspace=0)
    else:
        gs.update(wspace=0.08, hspace=0.25)

    matplotlib.rc('xtick', labelsize=6)
    matplotlib.rc('ytick', labelsize=6)
    fig.text(0.5, 0.01, 'Frequency f [1/d]', ha='center', fontsize=8)
    fig.text(0.05, 0.5, gls_obs.label["ylabel"], va='center',
             rotation='vertical', fontsize=8)

    if params['use_rms_as_jitter']:
        params['jitter'] = gls_obs.rms
    else:
        params['jitter'] = params['jitter']


    ylim_max = np.nanmax(gls_obs.power)
    axes = []
    save_var = np.empty((len(sim_freqs), params['mc_samples'], len(gls_obs.power)))
    for sim_freq_idx, freq in enumerate(sim_freqs):
        print(f'\n Simulating Freq. # {sim_freq_idx+1:.0f}\n')
        gls_sim_powers = np.zeros((params['mc_samples'], len(gls_obs.power)))
        freqs_array = np.zeros((params['mc_samples'], len(gls_obs.freq)))
        sim_phases = np.zeros((params['mc_samples'], len(peaks_data[2])))
        if params['num_processes'] > 1:
            pbar = tqdm(total=params['mc_samples'])
            def update(*a):
                pbar.update()
            with Pool(params['num_processes']) as p:
                res_temp = [p.apply_async(af_calc.sample_gls,
                                          args=(times, gls_obs, freq,  params['jitter'],
                                                rvs_err, params['fbeg'],  params['fend'],
                                                params['object_name'], peaks_data,
                                                params['search_phase_range']),
                                          callback=update)
                            for i in range(pbar.total)]
                for idx, result in enumerate(res_temp):
                    gls_sim_powers[idx], \
                     freqs_array[idx], \
                     sim_phases[idx] = result.get()
        else:
            for i in tqdm(range(0, params['mc_samples'])):
                gls_sim_powers[i], \
                 freqs_array[i], \
                 sim_phases[i] = af_calc.sample_gls(times, gls_obs, freq,
                                                    params['jitter'], rvs_err, params['fbeg'],
                                                    params['fend'], params['object_name'],
                                                    peaks_data,
                                                    params['search_phase_range'])
        ylim_max = np.nanmax([ylim_max, np.nanmax(np.percentile(gls_sim_powers,
                                                                97.5, axis=0))])
        if params['calc_metric']:
            metric=af_calc.get_metric(gls_obs=gls_obs,
                gls_sim_powers=gls_sim_powers)
            print("\n")
            print(f'rms={metric:.2f}')

        ax = af_plots.plot_panel_row(fig, gs, plot_row=sim_freq_idx,
                                     panel_borders=panels, sim_freq=freq,
                                     gls_obs=gls_obs,
                                     gls_sim_powers=gls_sim_powers,
                                     peak_pos=peaks_data[0],
                                     obs_phases=peaks_data[2],
                                     sim_phases=sim_phases,
                                     hide_xlabel=params['hide_xlabel'])
        axes.append(ax)
    #     save_var[idx] = gls_sim_powers
    # Table(save_var).write(os.path.join(params['savepath'], 'simulated_GLS.fits'),
    #                       overwrite=True)

    if params['plot_additional_period_axis']:
        fig.text(0.5, 0.96, 'Period [d]', ha='center', fontsize=8)
        for panel, base_panel in zip(axes[0], axes[-1]):
            ax2 = panel.twiny()
            ax2.set_xticks(base_panel.get_xticks())
            ax2.set_xticklabels(af_plots.convert_frequency_to_period(base_panel.get_xticks()))
            ax2.set_xlim(panel.get_xlim())
            ax2.spines['right'].set_visible(False)
            ax2.spines['left'].set_visible(False)
            ax2.set_yticks([])
            ax2.tick_params(labelright='off')
    for row in axes:
        for idx, ax in enumerate(row):
            ax.set_ylim(0, ylim_max*1.45)
            if idx == 0:
                ax.yaxis.set_major_locator(MaxNLocator(nbins=3,
                                                       prune='both',
                                                       min_n_ticks=3))
                ax.yaxis.set_label_position("left")
                ax.yaxis.tick_left()
            else:
                ax.set_yticks([])
                ax.tick_params(labelleft='off')

    plt.tight_layout()
    save_string = os.path.join(params['savepath'], params['object_name']
                               + '_'+str(params['test_period']).replace('.', 'p')+'d_'
                               + str(params['sampling_freq']).replace('.', 'p')
                               + 'd_alias_test.pdf')
    plt.savefig(save_string,
                bbox_inches='tight', dpi=600)
    plt.show()
    print(f'\n \n Plot saved to {save_string}')


if __name__ == "__main__":
    main()
