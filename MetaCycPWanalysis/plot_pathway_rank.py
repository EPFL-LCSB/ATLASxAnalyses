__author__ = 'anastasia'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

green = '#5a9957'
blue = '#1c74ac'
violet = '#efa6b4'

def plotPwRank(df, ax, color):
    df_scope = df.query('{}==1&{}==1&pw_rank<=5&pw_rank>0'.format(column_network, column_path))
    values = list(set(df_scope['pw_rank'].to_list()))
    counts = [df_scope['pw_rank'].to_list().count(i) for i in values]
    ranges = [(6, 10), (10, 16), (16, 21), (21, 51), (51, 101), (101, 301), (301, 501), (501, 1001), (1001, 5001), (5001, 10001)]
    for index, range in enumerate(ranges):
        low = range[0]
        high = range[1]
        df_scope = df.query('{}==1&{}==1&pw_rank<=@high&pw_rank>@low'.format(column_network, column_path))
        values.append(6+index)
        counts.append(len(df_scope))
    values.append(7+index)
    df_scope = df.query('{}==1&{}==1&pw_rank==0'.format(column_network, column_path))
    counts.append(len(df_scope))
    ax.bar(values, counts, color = color)
    t = np.arange(1, 17).tolist()
    ax.set_xticks(t)
    plt.yticks(np.arange(0, max(counts)+1, 50))
    t_labels = np.arange(1, 6).tolist()
    t_labels.extend(['6-10', '10-15', '16-20', '21-50', '51-100', '101-300', '301-500', '501-1000', '1001-5000', '5001-10000', '10000+'])
    ax.set_xticklabels(t_labels)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_ylabel('Number of pathways', fontsize=16)
    ax.set_xlabel('Rank of the exactly reconstructed MetaCyc pathway')
    ax.tick_params(axis='x', labelrotation = 90)

db_scopes = ['MetaCyc_hp','MetaCyc', 'chemATLAS_hp', 'chemATLAS']
titles = ['MetaCyc reactions', 'BNICE.ch-curated MetaCyc','All ATLASx', 'BNICE.ch-curated ATLASx']

for index, db_scope in enumerate(db_scopes):
    fig, ax = plt.subplots()

    column_network = db_scope+"_coveredByNetwork"
    column_path = db_scope+"_hasPath"
    metaCycPwDataset = pd.read_csv('output/pw_ranking_'+db_scope+'.csv')
    print(db_scope)

    ax.set_title(titles[index])
    plotPwRank(metaCycPwDataset, ax, green)


    fig.tight_layout()
    fig.set_size_inches(7, 4)
    plt.subplots_adjust(left=None, bottom=0.45, right=None, top=None, wspace=None, hspace=0.8)

    plt.savefig('plots/pathways_rank_{}.png'.format(db_scope))