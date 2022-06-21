#%% INFO
# This script plots a graph with BER curves from generated .txt files
# obtained for MC-CDMA having different SF cases,  for 2 taps.

#%% REQUIRED MODULES
# paths and regexp
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
#
fntsize = 14
plt.rc('font', size=fntsize)
plt.rc('legend', fontsize=fntsize-2)

#%% PARAMS
# save figures or not (show them only)
save = False

#%% LOAD DATA
BER = {}

data = np.loadtxt('BER_NoOffset_QPSK_MMSE_L1.txt')
BER_SF1 = np.array(data.T[1])
BER_SF1[BER_SF1 == 0] = np.nan
EbNo = data.T[0]

data = np.loadtxt('BER_TimingOffset_QPSK_MMSE_L1.txt')
BER_SF2 = np.array(data.T[1])
BER_SF2[BER_SF2 == 0] = np.nan

data = np.loadtxt('BER_FreqOffset_QPSK_MMSE_L1.txt')
BER_SF4 = np.array(data.T[1])
BER_SF4[BER_SF4 == 0] = np.nan


#%%
plt.figure()
plt.semilogy(EbNo, BER_SF1, label=r'$No Offset$')
plt.semilogy(EbNo, BER_SF2, label=r'$Timing Offset$')
plt.semilogy(EbNo, BER_SF4, label=r'$Frequency Offset$')

plt.grid()
#
plt.xlim([0, 40])
plt.ylim([1e-5, 1])
#
plt.xlabel(r'$E_b/N_0$ [dB]')
plt.ylabel(r'BER')
#
plt.legend()

# save or show figure
if save:
    fname= 'BER_Taps1_imperfection.png'
    plt.savefig(fname=fname, bbox_inches='tight')
    plt.close()
else:
    plt.tight_layout()
    plt.show()
