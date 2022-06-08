#%% REQUIRED MODULES
# paths and regexp
from math import log2
import os
import re
#
from scipy import special as sfun
import numpy as np
from numpy import pi
#
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D # for custom legend
#
# plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# plt.rc('text.latex', preamble=r'\usepackage{amsmath,amssymb}')
# matplotlib.verbose.level = 'debug-annoying'
#
fntsize = 14
plt.rc('font', size=fntsize)
plt.rc('legend', fontsize=fntsize-2)

#%% PARAMS
# fixed set of parameters for which (available) BER curves will be plotted
prms = r'BER_\w+_MMSE_L8.txt'

# save figures or not (show them only)
save = False

#%% LOAD DATA for plain OFDM
# get all file names in current folder
fnames = ' '.join(os.listdir('./'))

# extract .txt files with prms in their name (REGEXP)
fnames = re.findall(prms, fnames)

# modulations for which the results are available
mods_raw = re.findall(r'BER_(\w+)_MMSE_L8.txt', ' '.join(fnames))

# all potential modulations
mods_all = ['QPSK', '16QAM', '64QAM', '256QAM']

# corresponding bits/symbol
ldM = [2, 4, 6, 8]

# sort rates in an increasing order
mods = []
for mod in mods_all:
    if mod in mods_raw: mods.append(mod)

# lists to store data
snr_dB = []
BER = []

# load files
for mod in mods:
    fname = 'BER_{}_MMSE_L8.txt'.format(mod)
    data = np.loadtxt(fname)
    snr_dB.append(data.T[0])
    BER.append(data.T[1])

#%% LOAD DATA for OFDM with RATE ADAPTATION
num_taps = 8

data = np.loadtxt('BER_RA_MMSE_L{}.txt'.format(num_taps))
snr_ra_dB = data.T[0]
BER_ra = data.T[1]

data = np.loadtxt('RATE_RA_MMSE_L{}.txt'.format(num_taps))
RATE_ra = data.T[1]

#%% SNR in AWGN
# def bep_mqam_bnd(snr, M):
#     # Goldsmith (looser bound)
#     return 0.2*np.exp(-1.5*snr/(M-1))

# def qfun(x):
#     return  0.5*sfun.erfc(x/np.sqrt(2))

# def bep_mqam_exact(snr, M):
#     return (1 - (1 - 2*(1-1/np.sqrt(M)) * qfun(np.sqrt(3/(M-1)*snr)) )**2 )/np.log2(M)


# # theoretical BER in AWGN
# snr1_dB = np.linspace(EbNo_wf[0],EbNo_wf[-1],num=100)
# snr1 = 10**(snr1_dB/10)
# Pb_awgn = []
# for k in range(len(mods)):
#     Pb_awgn.append(bep_mqam_bnd(snr1, 2**ldM[k]))
#     # Pb_awgn.append(bep_mqam_exact(snr1, 2**ldM[k]))

#%% PLOT BER
Pb_th = 1e-3

plt.figure()
plt.axhline(Pb_th, ls=':', color='gray') # target BER
# BER for fixed modulations
for k in range(len(mods)):
    plt.semilogy(snr_dB[k], BER[k], '--', label=str(mods[k]))
    # SNR = EbNo[k] + 10*np.log10(ldM[k]) + 10*np.log10(cploss) #OLD
    # plt.semilogy(SNR, BER[k], '--', label=str(mods[k])) # vs SNR
# theoretical BER in AWGN
# plt.gca().set_prop_cycle(None) # reset color cycle
# for k in range(len(mods)):
#     plt.semilogy(snr1_dB, Pb_awgn[k], ':', label=str(mods[k]))
#
plt.semilogy(snr_ra_dB, BER_ra, 'k-', label='adaptive')
plt.grid()
#
plt.legend()
#
plt.xlim([0, 50])
plt.ylim([1e-5, 0.5])
#
plt.xlabel(r'SNR [dB]') # r'$E_b/N_0$ [dB]'
plt.ylabel(r'BER')
#
# save or show figure
if save:
    fname= 'BER_wf_ofdm.png'
    plt.savefig(fname=fname, bbox_inches='tight')
    plt.close()
else:
    plt.tight_layout()
    plt.show()

# %% PLOT THROUGHPUT
N_crr = 64; # number of carriers
plt.figure()
plt.gca().set_axisbelow(True)
for k in range(len(mods_all)):
    plt.plot(snr_ra_dB, ldM[k]*N_crr*np.ones_like(snr_ra_dB), '--', label=str(mods_all[k]))
plt.plot(snr_ra_dB, RATE_ra, 'k-', label='adaptive')
plt.grid()
#
plt.legend(loc='lower right')
#
plt.xlim([0, 40])
plt.ylim([0, np.log2(256)*N_crr])
#
plt.yticks([N_crr*x for x in ldM])
#
plt.xlabel(r'SNR [dB]') # r'$E_b/N_0$ [dB]'
# plt.ylabel(r'$\mathbb{E}(\mathrm{ld}M)$ [bits/symbol]')
plt.ylabel(r'Throughput [bits/block]')
#
# save or show figure
if save:
    fname= 'RATE_wf_ofdm.png'
    plt.savefig(fname=fname, bbox_inches='tight')
    plt.close()
else:
    plt.tight_layout()
    plt.show()

#%%