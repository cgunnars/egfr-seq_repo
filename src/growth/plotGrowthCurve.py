import argparse
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', help = 'clean data to plot')
parser.add_argument('-out', help = 'output file')
parser.add_argument('--axenic', action='store_true')
args = parser.parse_args()

data = pd.read_excel(args.i)

data = data[(data['axenic'] == args.axenic) & (data['method'] == 'lux') & (data['dose'] == 10)]
order = data[(data['day'] == data['day'].max())].groupby('condition').agg({'norm mfi': 'mean'}).sort_values(by=['norm mfi']).index
sns.set_style('white')
sns.set_palette('deep')

g = sns.FacetGrid(data, col = 'condition', hue = 'donor', col_order=order, col_wrap = 5)
g = g.map(sns.lineplot, 'day', 'norm mfi')
g = g.map(plt.axhline, ls=':', y=0.5, color='grey')
g = g.map(sns.scatterplot, 'day', 'norm mfi')

g.savefig(args.out)
