

##%matplotlib inline

import numpy as np
from scipy.io import loadmat, savemat

from fooof import FOOOFGroup, FOOOF

####from fooof.plts.periodic import plot_peak_fits, plot_peak_params
####from fooof.plts.aperiodic import plot_aperiodic_params, plot_aperiodic_fits
####
####from fooof.bands import Bands
####from fooof.analysis import get_band_peak_fg
####from fooof.utils import trim_spectrum
####from fooof.utils.data import subsample_spectra
####from fooof.sim.gen import gen_aperiodic
####from fooof.data import FOOOFSettings
####from fooof.plts.templates import plot_hist
####from fooof.plts.spectra import plot_spectra
####
##
### Load the mat file 
data = loadmat('power_spectra.mat')
##
### Unpack data from dictionary, and squeeze numpy arrays
freqs = np.squeeze(data['freqs']).astype('float')
psds = np.squeeze(data['avgpsd']).astype('float')
### ^Note: this also explicitly enforces type as float (type casts to float64, instead of float32)
###  This is not strictly necessary for fitting, but is for saving out as json from FOOOF, if you want to do that
##
### Transpose power spectra, to have the expected orientation for FOOOF
psds = psds.T
##
##
####
##### Initialize FOOOFGroup object
####fg = FOOOFGroup()
####fg.report(freqs, psds, [.2, 45])


fm = FOOOF()
fm.report(freqs, psds, [.2, 45])
fooof_results = fm.get_results()
fooof_results_dict = fooof_results._asdict()
savemat('fooof_results.mat', fooof_results_dict)

##fm.save('fit_data', save_data=True)
##exps = fm.get_params('aperiodic_params', 'exponent')
##savemat('exps.mat', {'exps' : exps})
##
####
##### Save out a specific FOOOF measure of interest - for example, slopes
####exps = fg.get_params('aperiodic_params', 'exponent')
####savemat('exps.mat', {'exps' : exps})
######
####### Save out fooof results to json file
#######  There is a utility file to load this json file directly into Matlab
####fg.save('fooof_results', save_results=True)
######
####
####fg.save('fit_data', save_data=True)
##
##
##
##
####
####fm = FOOOFGroup()
######fm.report(freqs, psds, [1, 30])
####fm.fit(freqs, psds, [3, 40])
####
####init_ap_fit = gen_aperiodic(fm.freqs, fm.get_params('aperiodic_params'))
##
