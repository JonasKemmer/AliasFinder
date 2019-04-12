import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Cursor
from af_utils import *
from detect_peaks import detect_peaks
from scipy.stats import circmean, circstd
import warnings
from tqdm import tqdm
warnings.filterwarnings("ignore")

verbose = False
object_name = 'GJ411'
rv_files = ['/home/jkemmer/Documents/CARMENES/projects/Gl411/Data/J11033+359.avc.dat']
test_freq = 12.940  # [days]

mc_samples = 10

panel_width = 1  # [days] width of the plot panels around the simulated frequencies
fbeg = 0.0001  # GLS frequency range (change to get better zoom, default is good for first overview
fend = 2.5
power_threshold = 0.06 # power threshold (ZK normalization) to find peaks in the periodograms
search_phase_range = 0.00005 # frequency range to search for new maxima in simulated data bases on real data GLS peaks

offsets = np.array([0])
substract = False  # remove sinusoidal singnal using sinus function 1: Subtract most significant period, 2 : Subtract second m-s-p
jitter = 0  # jitter of planet fit to insert in simulated data
use_rms_as_jitter = True  # Use rms as jitter (set to 1)


################################################################################
########################### START of PROGRAM ###################################
################################################################################
# read in rv data
times, rvs, rvs_err = read_rvs(rv_files, offsets)

# First analyze the spectral window function of the data and select a sampling
# frequency who's aliases should be analyzed
window_freqs = np.linspace(fbeg, fend, len(times)*30)
window_powers = get_spectral_window_function(window_freqs, times)

print('Please select the sampling frequency you want to test from the plot')
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.title('Spectral window function of the data\n'
          'Please select the sampling frequency you want to test '
          '(Press Spacebar+click)')
plt.xlabel('Frequency f [1/d]')
plt.ylabel(r'Power W($\nu$)')
ax.plot(window_freqs, window_powers)
plt.savefig(f'{object_name}_spectral_window_function.pdf',
            bbox_inches='tight', dpi=400)
cursor = Cursor(ax, useblit=True, color='k', linewidth=1)
zoom_ok = False
print('\nZoom or pan to view, \npress spacebar when ready to click:\n')
while not zoom_ok:
    zoom_ok = plt.waitforbuttonpress()
sampling_freq = plt.ginput(1)
sampling_freq = np.round(sampling_freq[0][0], 5)
plt.close()
print(f'Your selected sampling frequency is {sampling_freq}')


# Second, analyze the RV data
gls_obs = get_gls(times, rvs, rvs_err, fbeg, fend, object_name)

# if kepler singnal has to be substracted
if substract == 1:
    print('You chose to substract the most significant Period!')
    rv_array=rvs-gls_obs.sinmod(times)
    gls_obs = get_gls(times, rv_array, rvs_err, fbeg, fend, object_name,
                      verbose=verbose)
    rvs=rv_array
    print('Above Signal now significant!')

if substract == 2:
    print('You chose to substract the 2 most significant period!')
    rv_array = rvs - gls_obs.sinmod(times)
    gls_obs = get_gls(times, rv_array, rvs_err, fbeg, fend, object_name,
                      verbose=verbose)
    rv_array = rv_array-gls_obs.sinmod(times)
    gls_obs = get_gls(times, rv_array, rvs_err, fbeg, fend, object_name,
                      verbose=verbose)
    rvs=rv_array
    print('Above Signal now most significant!')

peaks_data = get_phase_info(gls_obs, power_threshold=power_threshold)
np.savetxt('{}_data.txt'.format(object_name), np.transpose(peaks_data),
           fmt='%s', newline='\n', delimiter='\t')

if use_rms_as_jitter:
    jitter  = gls_obs.rms

# Select the Periods which are to be simulated
sim_freqs = np.zeros(3)
for idx, m in enumerate([0, -1, 1]):
    while True:
        wanted_freq = np.round(calc_alias_freq(test_freq, sampling_freq, m), 5)
        closest_freq = np.round(find_nearest(peaks_data[0], wanted_freq), 5)
        check = input(f'Please check {idx+1}. frequency you want to simulate\n'
                      f'The predicted freq is: {wanted_freq}\n'
                      f'The closest found freq is:            {closest_freq}\n'
                      f'Do they fit? (yes/no)\n'
                      f'(If you want to manually choose another '
                      f'frequency please also type \'no\')  ')
        if check.lower() == 'yes' or check.lower() == 'y':
            sim_freqs[idx] = closest_freq
            break
        elif check.lower() == 'no' or check.lower() == 'n':
            print('The alias was not found in the data. Try to lower '
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
                    print('Your chosen frequency is not in the peak array. '
                          'Try again.')
            else:
                raise Exception('No frequency to simulate!')


# Definition of the borders of the plot panels
f1s = 1/(test_freq+panel_width)
f1e = 1/(test_freq-panel_width)
f2s = np.abs(1/(test_freq-panel_width)-sampling_freq)
f2e = np.abs(1/(test_freq+panel_width)-sampling_freq)
f3s = np.abs(1/(test_freq+panel_width)+sampling_freq)
f3e = np.abs(1/(test_freq-panel_width)+sampling_freq)
panels = np.array([[f1s, f1e], [f2s, f2e], [f3s, f3e]])#, [f2s+1, f2e+1], [f3s+1, f3e+1]])

# Initializing the plot and starting the simulations
fig = plt.figure(figsize=(7, 5))
gs = gridspec.GridSpec(3, 3)
gs.update(wspace=0.08, hspace=0.25)
matplotlib.rc('xtick', labelsize=9)
matplotlib.rc('ytick', labelsize=9)
fig.text(0.5, 0.04, 'Frequency f [1/d]', ha='center')
fig.text(0.04, 0.5, gls_obs.label["ylabel"], va='center', rotation='vertical')

for sim_freq_idx, freq in enumerate(sim_freqs):
    print(f' Simulating Freq. # {sim_freq_idx+1:.0f}\n')
    gls_sim_powers = np.zeros((mc_samples, len(gls_obs.power)))
    freqs_array = np.zeros((mc_samples, len(gls_obs.freq)))
    sim_phases = np.zeros((mc_samples, len(peaks_data[2])))
    for i in tqdm(range(0, mc_samples)):
        rvs_sim = np.ones(np.shape(times)) \
                   + sim_any_sinmode(gls_obs, freq, times) \
                   + np.random.normal(0, np.sqrt(jitter**2), times.size)
        ls_sim = get_gls(times, rvs_sim, rvs_err, fbeg, fend, object_name,
                         freq_array=gls_obs.freq)
        dummy_freq_array = np.zeros(np.size(peaks_data[0]))
        #search for phases of max power using a certian frequency range and the prior of data peaks
        for j in range(0, np.size(peaks_data[0])):
            index_frequencies = np.where(np.logical_and(ls_sim.freq >= peaks_data[0][j]
                                                         - search_phase_range,
                                                        ls_sim.freq <= peaks_data[0][j]
                                                         + search_phase_range))
            index_maxfreqs = max(np.arange(len(ls_sim.p[index_frequencies])),
                                 key=ls_sim.p[index_frequencies].__getitem__)
            index_maxfreq = np.argwhere(ls_sim.freq==ls_sim.freq[index_frequencies][index_maxfreqs])[0]

            dummy_freq_array[j] = ls_sim.freq[index_maxfreq]
        peaks_sim = get_phase_info(ls_sim, frequency_array=dummy_freq_array,
                                   sim=False)
        gls_sim_powers[i] = ls_sim.power
        freqs_array[i] = ls_sim.freq
        sim_phases[i] = (peaks_sim[2] % 1) * 2. * np.pi

    plot_panel_row(fig, gs, plot_row=sim_freq_idx, panel_borders=panels,
                   sim_freq=freq,
                   gls_obs=gls_obs, gls_sim_powers=gls_sim_powers,
                   peak_pos=peaks_data[0], obs_phases=peaks_data[2],
                   sim_phases=sim_phases)

plt.tight_layout()
plt.savefig('GJ411_alias_test.pdf', bbox_inches='tight', dpi=400)
plt.show()
