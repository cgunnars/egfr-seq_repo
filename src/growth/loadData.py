import os
import pandas as pd
import numpy as np
import argparse
import re


# conditions in format cond1 cond2 
# ex. sm drug // feso4 ima

## TODO: deal with uninf, dmso, media
## takes in row['value_x'] = sm drug
##          row['value_y'] = feso4 ima
## sets the column names to be sm: feso4, drug: ima
def set_cond(row, condnames, ctrls=['uninf','dmso','media']):
    colnames = row['condition_name'].split(' ')
    colvals = row['condition'].split(' ')    
    
    for colname, colval in zip(colnames, colvals):
       if colname in condnames:
            row[colname] = colval

    #remove donor from condition to avoid issues later  
    row['condition'] = ' '.join([x for x, y in zip(colvals, colnames) if y != 'donor'])
    return row

def walkdir(dir, nrow = 6):
    #preallocate mem
    data = pd.DataFrame()
    plate = pd.DataFrame()
    condnames = []

    ignore_names = ['~', '.DS_Store', '.git', '.ipynb', '-data.xlsx', '.txt'] 
    for i, name in enumerate(os.listdir(dir)):
        path = os.path.join(dir, name)
        if (not any(s in path for s in ignore_names)):
            if os.path.isdir(path):
                data = data.append(walkdir(path), ignore_index = True)
            
            elif ('PlateLayout' in path):
                temp = pd.read_excel(path, header = None, index_col = None)
                temp = temp.dropna(axis=0, how='all').dropna(axis=1, how='all')
                #temp.columns = [int(col) + 1 for col in temp.columns]
                # Plate setup:
                # Drug info w/ spaces separating sm/drug/dose
                # Metadata (is this a sm/drug/dose)
                # What to normalize to
                temp_condition = reindex(temp.iloc[0:nrow, :], nrow = nrow).drop(columns=['variable'])
                temp_condition = temp_condition.rename(columns = {'value': 'condition'})

                temp_metadata = reindex(temp.iloc[nrow:nrow*2, :], nrow = nrow).drop(columns=['variable'])
                temp_normto = reindex(temp.iloc[nrow*2:nrow*3, :], nrow = nrow).drop(columns=['variable'])
                
                # by default, say the condition is a drug
                # otherwise, need to rename the condition according the metadata
                if temp_metadata.empty:
                    temp = temp_condition
                    temp.rename(columns={'condition': 'drug'})
                else:
                    temp_metadata = temp_metadata.rename(columns = {'value': 'condition_name'})
                    temp = pd.merge(temp_condition, temp_metadata, on='position', how='outer')
                    temp_condnames = temp_metadata['condition_name'].str.split(expand=True).values.flatten()
                    temp_condnames = np.unique([x for x in temp_condnames if (x is not None and x != 'XXX')])
                    temp = temp.apply(lambda x: set_cond(x, temp_condnames), axis=1)
                    temp = temp.drop(columns=['condition_name'])
                
                # by default, normalize to dmso
                # otherwise, use metadata 
                if temp_normto.empty:
                    temp['norm_to'] = 'dmso'
                else:
                    temp_normto = temp_normto.rename(columns = {'value': 'norm_to'})
                    temp['norm_to'] = temp_normto['norm_to']
                # get names of subconditions
                temp['platename'] = name
                plate = plate.append(temp, ignore_index = True)
            else:
                temp = getdata(path)
                temp['day'] = getdate(path)
                temp['method'] = getmethod(path)
                temp['donor'] = getdonor(path)
                temp['exp'] = dir
                temp['axenic'] = ('axenic' in path)
                temp['filename'] = name
                data = data.append(temp, ignore_index = True)
    if 'variable' in data.columns: 
        data = data.drop(columns=['variable'])
        
    if (not plate.empty):
        if (len(plate['platename'].unique()) > 1): 
            temp = pd.DataFrame()
            #slice into matching plate / layout combos, merge, and accumulate
            for platename in plate['platename'].unique():
                platename = platename.split('_')[1]
                
                temp = temp.append(pd.merge(
                                        data[data['filename'].str.contains(platename)], 
                                        plate[plate['platename'].str.contains(platename)], 
                                        on='position'),
                                   ignore_index = True)
                
            data = temp
            print(platename)
        else:
            data = pd.merge(data, plate, on='position')
    return data

def getdata(path):
    if ('CJ_LDH' in path):
        skiprows = 62
        data = reindex(pd.read_excel(path, skiprows=skiprows, header=None, usecols = list(np.arange(1, 13)), nrows = 6)) 
    elif ('BB_alamar' in path):
        data_570 = reindex(pd.read_excel(path, skiprows=35, header=None, usecols = list(np.arange(1,13)), nrows=6))
        data_570 = data_570.rename(columns={'value': 'abs570'}).drop(columns='variable')
        data_600 = reindex(pd.read_excel(path, skiprows=59, header=None, usecols = list(np.arange(1,13)), nrows=6))
        data_600 = data_600.rename(columns={'value': 'abs600'}).drop(columns='variable')
        data = pd.merge(data_570, data_600, on='position')
    else:
        skiprows = 35
        df = pd.read_excel(path, skiprows=skiprows, header=None, usecols = list(np.arange(1, 13)), nrows = 6)
        data = reindex(df)
    return data

def reindex(data, nrow = 6):
    data.index = data.index + 2
    data.columns = [chr(int(i) + 65) for i in data.columns]
    data = data.melt()
    data['position'] = [i + str(j + 2) for i, j in zip(data['variable'], np.remainder(data.index.values, nrow))]
    return data

## works for /XXX_d3_ /XXX_D3_ /XXX_day3_ /XXX_Day3_ or at end of file name (e.g. /XXX_d3.xlsx)
def getdate(path):
    name = path.split(os.path.sep)[-1]
    date = re.search(r'((?<=_[dD])|(?<=_[dD]ay))\d(?=[_.])', name)
    if date:
        date = date[0]
    return date

def getdonor(path):
    name = path.split(os.path.sep)[-1]
    donor = re.search(r'((?<=_[dD]onor)|(?<=_))[A-Za-z](?=[_.])', name)
    if donor:
        donor = donor[0]
    return donor

def getdose(path):
    name = path.split(os.path.sep)[-1]
    dose = re.search(r'(?<=dose)(.*?)(?=[_.])', name)
    if dose:
        dose = float(re.sub('pt', '.', dose[0]))
    else:
        dose = 10
    return dose

def getmethod(path):
    methods = ['LDH', 'lux', 'OD600', 'alamar']
    name = path.split(os.path.sep)[-1]
    method = name.split('_')[1]
    if method not in methods:
        method = 'lux' #CJ_LDH_, GHB_lux_, AKB_OD600
    return method

parser = argparse.ArgumentParser()
parser.add_argument('-dir', help = 'name of main directory to traverse')
parser.add_argument('-out', help = 'name of output file')
args = parser.parse_args()

data = walkdir(args.dir)
data = data.rename(columns={'value': 'raw mfi'})

if (any(data.columns.str.contains('donor_x')) and not data['donor_x'].any()):
    data = data.rename(columns={'donor_y': 'donor'})
    print('dropped')
data.to_excel(os.path.join(args.dir, args.out))

