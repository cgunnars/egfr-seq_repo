import cirpy
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit import DataStructs
import pandas 
import seaborn
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import sklearn.metrics
import sklearn.linear_model

pandas.set_option('display.max_rows', 500)
pandas.set_option('display.max_columns', 30)
seaborn.set_theme(rc={'font.size': 12, 'axes.labelsize': 14, 'legend.fontsize': 10})

# input: data table containing SMILES strings for EGFRis of interest
# output: plot - chemical structure diagram for all EGFRis, aligned by shared scaffold
#         data table - matrix of tanimoto similarity (returned and written to file)
def getTanimoto(chem_info):
    cas = chem_info.loc[:, 'CAS#'] #['231277-92-2', '439081-18-2', '383432-38-0', '267243-28-7']
    smiles = chem_info.loc[:, 'SMILES'].values#[cirpy.resolve(x, 'smiles') for x in cas]
    names = chem_info.loc[:, 'Drug'].values


    ms = [Chem.MolFromSmiles(x) for x in smiles]
    for m, name in zip(ms, names):
        m.SetProp('_Name', name)


    q_smiles = 'C1=CC=C2C(=C1)C=NC=N2'
    q_core = Chem.MolFromSmiles(q_smiles)
    q_core = Chem.RemoveHs(q_core)
    subms = [x for x in ms if x.HasSubstructMatch(q_core)]
    highlight_q = [x.GetSubstructMatch(q_core) for x in ms]
    AllChem.Compute2DCoords(q_core)

    for m in subms:
        _ = AllChem.GenerateDepictionMatching2DStructure(m, q_core)

    opts = Draw.MolDrawOptions()
    opts.legendFraction = 0.1
    opts.minFontSize = 10
    opts.legendFontSize = 20
    opts.highlightRadius = 0.75 
    img=Draw.MolsToGridImage(ms, molsPerRow=5, subImgSize=(250, 250), highlightAtomLists=highlight_q, 
                             legends=[x.GetProp("_Name") for x in ms], drawOptions=opts)
    img.save('./fig/chem_info/scaffold_grid.pdf')


    fpgen = AllChem.GetRDKitFPGenerator()
    fps = [fpgen.GetFingerprint(x) for x in ms]

    qu, ta, sim = [], [], []

    # compare all fp pairwise without duplicates
    for n in range(len(fps)-1): # -1 so the last fp will not be used
        s = DataStructs.BulkTanimotoSimilarity(fps[n], fps[n+1:]) # +1 compare with the next to the last fp
        print(smiles[n], smiles[n+1:]) # which mol is compared with what group
    # collect the SMILES and values
        for m in range(len(s)):
            qu.append(names[n])
            ta.append(names[n+1:][m])
            sim.append(s[m])

    d = {'query':qu, 'target':ta, 'Similarity':sim}
    df_final = pandas.DataFrame(data=d)
    df_final = df_final.sort_values('Similarity', ascending=False)

    df_final.to_csv('./data/tanimoto.csv')

    tanimoto_matrix = df_final.pivot_table(index='query', columns='target', values='Similarity')
    tanimoto_matrix = tanimoto_matrix.combine_first(tanimoto_matrix.T).fillna(1)

    return(tanimoto_matrix)

# input: filename containing cleaned data
# output: median normalized control per condition at endpoint
def getCtrl(data='./data/lux_data/fig1/clean-data.xlsx'):

    ## read in control data
    ctrl_df      = pandas.read_excel(data)
    ctrl_df      = ctrl_df.loc[ctrl_df['day'] == 5, :]
    ctrl_summary = ctrl_df.groupby(['condition', 'donor']).agg({'norm mfi': np.mean})
    ctrl_summary = ctrl_summary.groupby(['condition']).agg({'norm mfi': np.median})

    return(ctrl_summary)

def calcCumulativeInhibition(data):
    ## set values with >100% activity as 100 to avoid going negative
    inhib = [100 - x if x < 100 else 0 for x in data]
    # gini curve y
    inhib_y = np.cumsum(np.sort(inhib) / np.sum(inhib))
    inhib_x = np.cumsum(np.ones(len(inhib)) / len(inhib))
    auc = sklearn.metrics.auc(inhib_x, inhib_y)

    gini = 1 - 2 * auc
    return(gini)
# returns colors for log2 control (e.g. -1 = 0.5, 0 = 1) ranging from green to pink
# ASSUMES YOU'RE ROUNDING TO NEAREST 0.5 of log2 normalized control
def getCtrlColors():
    ## hard code colors for now..... 
    bnd_colors = [-3.5, -2.5, -1.5, -1, -0.5, 0, 0.5, 1.5]
    rgb_colors = [(39 / 255, 100 / 255, 25 / 255), 
                  (77 / 255, 146 / 255, 33 / 255), 
                  (127 / 255, 188 / 255, 65 / 255), 
                  (127 / 255, 188 / 255, 65 / 255), 
                  (230 / 255, 245 / 255, 208 / 255), 
                  (1, 1, 1), 
                  (253 / 255, 224 / 255, 239 / 255), 
                  (222 / 255, 119 / 255, 174 / 255)]
    lut = dict(zip(bnd_colors, rgb_colors))
    return(lut)

def getInhib(chem_info, data='./data/chem_info/2011-nbt_inhibition.xls'):
    inhib_data = pandas.read_excel(data, header=1, index_col=0).transpose()
    inhib_data = inhib_data.set_index('compound CAS#:')
    
    # match on compound cas number
    ki_data = chem_info.set_index('CAS#')
    compounds = ki_data.index.dropna()
    compounds = (compounds[compounds.isin(inhib_data.index)])
    
    inhib_data = inhib_data.loc[compounds,:].fillna(np.nan)
    inhib_data.index = ki_data.loc[compounds,'Drug']
    return(inhib_data)

## read in chemical information 
chem_info = pandas.read_excel('./data/chem_info/EGFR_KIs-20220415.xlsx')
## MAKE CHEMICAL SIMILARITY PLOTS -- PUBCHEM
# plot clustered tanimoto distance, annotated according to compound intrabacterial effect
tm   = getTanimoto(chem_info)
ctrl = getCtrl()
ctrl_plot = ctrl.loc[tm.index.str.lower(), :]
lut  = getCtrlColors()
lognorm_ctrl = round(np.log2(ctrl_plot['norm mfi']) * 2) / 2
row_colors   = lognorm_ctrl.map(lut).values
col_colors   = row_colors

plt.figure(figsize=(3,3))
hm = seaborn.clustermap(tm, cmap='viridis', row_colors=row_colors, col_colors=col_colors)
hm.savefig('./fig/chem_info/tanimoto.svg', bbox_inches='tight')


### MAKE INHIBITION SPECIFICITY PLOTS -- ANASTASSIADIS
inhib_data = getInhib(chem_info)
# reorder control values to match plot order
lognorm_ctrl = lognorm_ctrl[inhib_data.index.str.lower()]
row_colors = lognorm_ctrl.map(lut).values
# subset to only compounds with strong inhibition
inhib_data_sparse = inhib_data.loc[:, (inhib_data[inhib_data < 75].count() > 1)]
plt.figure(figsize=(0.5*len(inhib_data_sparse.columns), 0.75 * len(inhib_data_sparse.index)))
seaborn.clustermap(data = inhib_data_sparse, mask = inhib_data_sparse.isnull(), cmap='viridis_r', col_cluster=False, row_cluster=False, row_colors=row_colors, cbar_pos=(0, .2, .03, .4))#, dendrogram_ratio=0)
plt.savefig('./fig/chem_info/EGFR-spec-heatmap.svg', bbox_inches='tight')



gini = inhib_data.apply(calcCumulativeInhibition, axis=1)


kd_data = pandas.read_excel('./data/chem_info/Klaeger_s3.xlsx', sheet_name=2, index_col=0).transpose()
compounds = [x for x in chem_info.loc[:,'Drug'] if x in kd_data.index]
kd_data = kd_data.loc[compounds, :]
kd_data.index = kd_data.index.str.lower()
kd_data = kd_data.join(ctrl)

plt.figure()
kdplot = seaborn.scatterplot(data=kd_data, x='norm mfi', y='EGFR')
plt.yscale('log')
kdplot.figure.savefig('./fig/chem_info/EGFR-kd-ctrl_a.pdf', bbox_inches='tight')



catds_data = pandas.read_excel('./data/chem_info/Klaeger_S6.xlsx', sheet_name=3, index_col=0, header=4)
catds_data = catds_data.loc[compounds, :]#['Gini', 'Selectivity', 'CATDS single']]
catds_data.index = catds_data.index.str.lower()

catds_data = catds_data.join(ctrl)


plt.figure(figsize=(3,3))
kdplot = seaborn.scatterplot(data = catds_data, x = 'norm mfi', y = 'Concentration [nM]', hue='CATDS single', size='CATDS single')
for line in range(0,catds_data.shape[0]):
     plt.text(
          catds_data["norm mfi"][line]+0.05,
          catds_data["Concentration [nM]"][line],
          catds_data.index[line],
          ha='left'
     ) ## label points with each drug
plt.yscale('log')
#plt.xscale('log', base=2)
x_lux = catds_data.loc[:, 'norm mfi'].values
y_kd  =  catds_data.loc[:, 'Concentration [nM]'].values


kdplot.figure.savefig('./fig/chem_info/EGFR-kd-ctrl.svg', bbox_inches='tight')
