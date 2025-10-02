import os 
import pandas as pd
import numpy as np
import argparse

# returns df with donor information
def assignDonor(data):
    data.loc[data['axenic'], 'donor'] = 'z'
    data['donor'] = data['donor'].fillna('a') #where no donor info assume donor A 
    data['donor'] = (data['donor'] + data['exp']).astype('category').cat.codes
    return data

# returns df with background (uninf for lux / carrier control for ldh) subtracted
def subBackground(data, methodname):
    if methodname == 'alamar': #data is grouped by donor
        O600 = 117216 
        O570 = 80586
        data['norm mfi'] = data['abs570']*O600 - data['abs600']*O570
    else:
        if methodname == 'ldh': #data is grouped by donor, then norm_to
            bgname = data.name[-1] #subtract absorbance of vehicle control 
        else: #assume lux
            bgname = 'uninf'
        bg = data[data['condition'].str.contains(bgname, na=False)]['raw mfi'].mean()
        if np.isnan(bg):
            bg = 0
        data['norm mfi'] = data['raw mfi'] - bg
    return data

# scale raw mfi to raw mfi of control at endpt (day 5)
def scale(data, methodname, endpt=5):
    if methodname == 'ldh':
        ctrlname = 'pc'
    else: #lux / alamar are grouped by normto
        ctrlname = data.name[-1]
    ctrl = data[data['day'] == min(data['day'].max(), endpt)]
    ctrl = ctrl[ctrl['condition'].str.contains(ctrlname, na=False)]
    if ctrl['norm mfi'] is not None:
        data['norm mfi'] = data['norm mfi'] / ctrl['norm mfi'].mean()
    return data

# normalize data on a per-donor basis by subtracting background then
# scaling according to endpt mfi of ctrl (given by df metadata)
# separate ldh data from all other data due to different normalization scheme
def normalize(data, group):
    data_ldh = data[data['method'] == 'ldh']
    data_lux = data[data['method'] == 'lux']
    data_alamar = data[data['method'] == 'alamar']
   
    group_dn = ['donor', 'norm_to']
    group_d  = ['donor']
    
    if not data_lux.empty:
        data_lux    = data_lux.groupby(group_d, group_keys=False, as_index=False).apply(lambda x: subBackground(x, 'lux')) 
    if not data_ldh.empty:
        data_ldh    = data_ldh.groupby(group_dn, group_keys=False, as_index=False).apply(lambda x: subBackground(x, 'ldh'))
    if not data_alamar.empty:
        data_alamar = data_alamar.groupby(group_d, group_keys=False, as_index=False).apply(lambda x: subBackground(x, 'alamar'))

    if group:
        group_dn = group_dn + [group]
        group_d  = group_d  + [group]

    data_lux    = data_lux.groupby(group_dn, group_keys=False, as_index=False).apply(lambda x: scale(x, 'lux'))
    data_ldh    = data_ldh.groupby(group_d,  group_keys=False, as_index=False).apply(lambda x: scale(x, 'ldh'))
    data_alamar = data_alamar.groupby(group_dn, group_keys=False, as_index=False).apply(lambda x: scale(x, 'alamar', 6))
        
    data = pd.concat([data_lux, data_ldh, data_alamar])
    return data

parser = argparse.ArgumentParser()
parser.add_argument('-dir', help = 'name of directory to process')
parser.add_argument('-i', help = 'name of raw input file')
parser.add_argument('-out', help = 'name of output file')
parser.add_argument('-group', help = 'extra group', required=False)
args = parser.parse_args()

data = pd.read_excel(os.path.join(args.dir, args.i))
data = data.applymap(lambda s: s.lower() if isinstance(s, str) else s) #
data = data[~(data['condition'].isin(['xxx', 'xx']))]
data = assignDonor(data)
data = normalize(data, args.group)
if('dose' in data.columns):
    data['dose'] = data['dose'].fillna(10) #if no dose info assume it is 1X
else:
    data['dose'] = 10
data.to_excel(os.path.join(args.dir, args.out))
