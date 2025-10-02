import argparse
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
sns.set_style('white')

parser = argparse.ArgumentParser()
parser.add_argument('-i', help = 'clean data to plot')
parser.add_argument('-out', help = 'output file')
parser.add_argument('-condlist', help = 'txt file containing list of conditions to plot', required = False)
parser.add_argument('-method', help='method to plot', required = False, default = 'lux')
parser.add_argument('-group', help = 'how to group data', required = False, nargs='+')
parser.add_argument('--lowthresh', help = 'plot line to show', required = False, action = 'store_false')
parser.add_argument('--highthresh', help = 'plot line to show', required = False, action = 'store_true')
parser.add_argument('--axenic', help = 'plot axenic data', required = False, action = 'store_true')
parser.add_argument('--tecreps', help = 'plot technical replicates', required = False, action = 'store_true')
parser.add_argument('-endpt', help='day of endpt', required = False, type=int, default = 5)
args = parser.parse_args()

ctrls = ['uninf', 'uninfected', 'media']
data = pd.read_excel(args.i)
endpt = data[(data['day'] == min(args.endpt, data['day'].max())) &
             (~data['condition'].isin(ctrls)) &
             (data['method'] == args.method)  &
             (data['axenic'] == args.axenic)]
condlist = None
if (args.condlist):
    with open(args.condlist, 'r') as fp:
        condlist = [line.strip() for line in fp]
    endpt = endpt[endpt['condition'].isin(condlist)]

cols = ['donor', 'condition']
if args.group:
    xvar = args.group[0]
    hue = args.group[1]
    striphue = args.group[1]
    cols = list(set(cols + args.group)) #remove duplicate groups
else:
    xvar='condition'
    hue = None
    striphue = None
#average technical replicates
if args.tecreps:
    endpt_tecmean = endpt
else:
    endpt_tecmean = endpt.groupby(cols, as_index=False).agg({'norm mfi': 'mean'})

endpt_allmean = endpt_tecmean.groupby(xvar).median()

print(endpt_allmean.head())
print(endpt_tecmean.head())
order = endpt_allmean.sort_values(by=['norm mfi'])

plt.figure(figsize=(6,3))
g = sns.boxplot(data=endpt_tecmean, x=xvar, y='norm mfi', hue = hue, 
                order = order.index, showfliers = False)
g = sns.stripplot(data=endpt_tecmean, x=xvar, y='norm mfi', hue = striphue, 
                  order = order.index, linewidth=1, dodge=True)
g.set_xticklabels(g.get_xticklabels(), rotation = 90)
if args.lowthresh:
    g.axhline(0.5, ls=':', color='grey')
if args.highthresh:
    g.axhline(1.50, ls=':', color='grey')
if args.group: #keep first legend for boxplot, delete duplicate stripplot legend
    ncond = max(endpt_tecmean.groupby(xvar).apply(lambda x: len(x[hue].unique())))
    handles, labels = g.get_legend_handles_labels()
    lgd = g.legend(handles[0:ncond], labels[0:ncond],
                   loc='upper left',
                   fontsize='large',
                   handletextpad=0.5)

g.spines['right'].set_visible(False)
g.spines['top'].set_visible(False)            
g.figure.savefig(args.out, bbox_inches='tight', width=6, height=2.5)

