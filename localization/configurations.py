import numpy as np

from pycbc.detector import add_detector_on_earth
from pycbc import psd
from simple_pe.detectors import Network

import logging
_logger = logging.getLogger('PESummary')
_logger.setLevel(logging.CRITICAL + 10)

df = 1./8
f_high = 4096
flen = int(f_high/df) + 1

flow = {'ET': 3,
        'CE': 5.2,
        'CE20': 5.2,
        'A+': 10,
        'A#': 10
        }

asharp = psd.from_txt('../ligo_strain/Asharp_strain.txt',
                      flen, df, flow['A#'])
ce40 = psd.from_txt('../ce_strain/cosmic_explorer_strain.txt',
                    flen, df, flow['CE'])
ce20 = psd.from_txt('../ce_strain/cosmic_explorer_20km_strain.txt',
                    flen, df, flow['CE'])
et = psd.analytical.EinsteinTelescopeP1600143(flen, df, flow['ET'])

add_detector_on_earth(name='C4', longitude=np.radians(-125), latitude=np.radians(46), 
                                    xangle=np.radians(260), yangle=np.radians(350) )
add_detector_on_earth(name='C2', longitude=np.radians(-94), latitude=np.radians(29), 
                                    xangle=np.radians(200), yangle=np.radians(290) )

psds = {
    'H1':asharp,
    'L1':asharp,
    'I1':asharp,
    'C4':ce40,
    'C2':ce20,
    'E1':et,
    'E2':et,
    'E3':et,
}

fmins = {'H1':flow['A#'],
    'L1':flow['A#'],
    'I1':flow['A#'],
    'C4':flow['CE'],
    'C2':flow['CE'],
    'E1':flow['ET'],
    'E2':flow['ET'],
    'E3':flow['ET'],
    
}

configurations = {}

configurations['HLA'] = ['H1', 'L1', 'I1']
configurations['40LA'] = ['C4', 'L1', 'I1']
configurations['40LET'] =  ['C4', 'L1', 'E1', 'E2', 'E3']
configurations['2040A'] =  ['C4', 'C2', 'I1']
configurations['2040ET'] =  ['C4', 'C2', 'E1', 'E2', 'E3']


def generate_networks(net_thresh=10, found_thresh=5, loc_thresh=4):
    networks = {}
    for key,ifos in configurations.items():
        networks[key] = Network(net_thresh)
        networks[key].generate_network_from_psds(ifos, psds, fmins, 
                                                 found_thresh=found_thresh, 
                                                 loc_thresh=loc_thresh)
    return networks
        