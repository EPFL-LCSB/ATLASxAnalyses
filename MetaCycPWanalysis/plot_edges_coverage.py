__author__ = 'anastasia'
__author__ = 'anastasia'
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

coverage_file = 'output/networkEdgesCoverage.csv'

green = '#5a9957'
blue = '#1c74ac'
violet = '#efa6b4'

# creating the directory for plots
if not os.path.exists('plots'):
    os.mkdir('plots')

def MetaCycVsATLASxAll():

    df = pd.read_csv(
        coverage_file)
    df = df.drop_duplicates()

    fig, ax = plt.subplots()

    bar1 =  ax.barh([0.02, 0.2],len(df),  0.1, color = violet)

    bar2, _ = plotNumberOfReconstructedPerNTDT(df, ['MetaCyc_hp_coveredByNetwork', 'chemATLAS_hp_coveredByNetwork'], ax, [green, green])

    labels = ['MetaCyc reactions', 'All ATLASx']
    ax.set_xlabel('Number of pathways', fontsize=14)
    ax.axes.yaxis.set_visible(False)
    ax.text(-50, 0.2,  labels[1], va='center', ha='center', fontsize=14, rotation = 90)
    ax.text(-50, 0.02, labels[0], va='center', ha='center', fontsize=14, rotation = 90)

    ax.text(len(df)+10, 0.1,  'Total: ' + str(len(df))+ ' pathways', ha='left', va='center', fontsize=14, rotation = 90)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlim(0, 1250)
    ax.set_ylim(-0.05, 0.27)
    fig.legend((bar1, bar2), ('Pathway not\nexactly reconstructed', 'Pathway exactly\nreconstructed'), loc=[0.10, 0.47], fontsize=14)
    plt.vlines(len(df), -0.03, 0.25, linestyles='solid')
    fig.tight_layout()
    fig.set_size_inches(4, 6)
    plt.gca().invert_yaxis()
    plt.savefig('plots/MetaCycVsATLASx.png', bbox = [6, 9])

def MetaCycVsATLASxBNICE():

    df = pd.read_csv(
        coverage_file)
    df = df.drop_duplicates()

    fig, ax = plt.subplots()

    bar1 =  ax.barh([0.02, 0.2],len(df),  0.1, color = violet)

    bar2, _ = plotNumberOfReconstructedPerNTDT(df, ['MetaCyc_coveredByNetwork', 'chemATLAS_coveredByNetwork'], ax, [green, green])

    labels = ['BNICE.ch-curated\nMetaCyc', 'BNICE.ch-curated\nATLASx']
    ax.set_xlabel('Number of pathways', fontsize=14)
    ax.axes.yaxis.set_visible(False)
    ax.text(-50, 0.2,  labels[1], va='center', ha='center', fontsize=14, rotation = 90)
    ax.text(-50, 0.02, labels[0], va='center', ha='center', fontsize=14, rotation = 90)

    ax.text(len(df)+10, 0.1,  'Total: ' + str(len(df))+ ' pathways', ha='left', va='center', fontsize=14, rotation = 90)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlim(0, 1250)
    ax.set_ylim(-0.05, 0.27)
    fig.legend((bar1, bar2), ('Pathway not\nexactly reconstructed', 'Pathway exactly\nreconstructed'), loc=[0.10, 0.45], fontsize=14)
    plt.vlines(len(df), -0.03, 0.25, linestyles='solid')
    fig.tight_layout()
    fig.set_size_inches(4, 6)
    plt.gca().invert_yaxis()
    plt.savefig('plots/MetaCycVsATLASxBNICE.png', bbox = [6, 9])

def plotNumberOfReconstructedPerNTDT(df, values, ax, color):

    counts = [df[i].to_list().count(1) for i in values]

    width = 0.1  # the width of the bars

    rects1 = ax.barh([0.02, 0.2],counts,  width, color = color)
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_xlabel('Number of pathways')

    autolabel(rects1, ax)

    return rects1

def autolabel(rects, ax):
        #Attach a text label above each bar in *rects*, displaying its height.
        for rect in rects:
            width = rect.get_width()
            print(width)
            ax.annotate('{}'.format(width),
                        xy=(width-110, rect.get_y() + rect.get_height() - 0.02),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=14)

MetaCycVsATLASxAll()
MetaCycVsATLASxBNICE()


######## Supplementary material plots #########
def plotLengthDistributionAll():
    df = pd.read_csv(coverage_file)
    df = df.drop_duplicates()
    df['Length'] = df.apply(calculate_length, axis=1)

    fig, ax = plt.subplots(figsize=(16,9))
    bar1 = plotLengthTotal(df,  ax, violet)
    bar4 = plotLengthOfReconstructedPerNTDT(df, ['chemATLAS_coveredByNetwork','chemATLAS_hp_coveredByNetwork'], ax, green) #['chemATLAS_edges','chemATLAS_hp_edges']

    ax.set_ylabel('Number of pathways', fontsize=20)
    ax.set_xlabel('Pathway length', fontsize=20)
    fig.legend((bar1, bar4), ('Pathway not\nexactly reconstructed', 'Pathway exactly\nreconstructed'), loc=[0.75, 0.75], fontsize=18)
    plt.xticks(np.arange(0, 30, 5), fontsize = 18)
    plt.yticks(np.arange(0, 450, 25), fontsize = 18)
    plt.xlim(0, 27)
    plt.ylim(0, 400)
    fig.tight_layout()
    plt.savefig('plots/plotMetacycPathwayLengthDistributionAll_all.png', bbox = [9, 9])

def plotLengthDistributionAllCrop():
    df = pd.read_csv(coverage_file)
    df = df.drop_duplicates()
    df['Length'] = df.apply(calculate_length, axis=1)

    fig, ax = plt.subplots(figsize=(9,3))
    _ = plotLengthTotal(df,  ax, violet)
    _ = plotLengthOfReconstructedPerNTDT(df, ['chemATLAS_coveredByNetwork','chemATLAS_hp_coveredByNetwork'], ax, green) #['chemATLAS_edges','chemATLAS_hp_edges']

    fig.tight_layout()
    plt.xticks(np.arange(10, 30, 5), fontsize = 21)
    plt.yticks(np.arange(5, 11, 5), fontsize = 21)
    plt.xlim(9.5, 27)
    plt.ylim(0, 12)
    plt.savefig('plots/plotMetacycPathwayLengthDistributionAll_crop.png', bbox = [3, 9])

def plotLengthOfReconstructedPerNTDT(df, search_setups, ax, color):

    setup1 = search_setups[0]
    setup2 = search_setups[1]
    values1 = list(set(df['Length'].to_list()))

    width = 0.4  # the width of the bars

    counts1 = [df[df[setup1]==1]['Length'].to_list().count(i) for i in values1]

    values1 = [i+ width / 2 for i in values1]
    rects1 = ax.bar(values1, counts1, width, label=setup1, color=color)

    values2 = list(set(df['Length'].to_list()))

    counts2 = [df[df[setup2]==1]['Length'].to_list().count(i) for i in values2]

    values2 = [i - width / 2 for i in values2]
    rects2 = ax.bar(values2, counts2, width, label=setup2, color=color)

    return rects1, rects2

def plotLengthTotal(df,  ax, color):

    values1 = list(set(df['Length'].to_list()))

    width = 0.4  # the width of the bars

    counts1 = [df['Length'].to_list().count(i) for i in values1]

    values1 = [i+ width / 2 for i in values1]
    rects1 = ax.bar(values1, counts1, width, color=color)

    values2 = list(set(df['Length'].to_list()))
    counts2 = [df['Length'].to_list().count(i) for i in values2]

    values2 = [i - width / 2 for i in values2]
    rects2 = ax.bar(values2, counts2, width, color=color)

    return rects1, rects2

def calculate_length(row):
    return len(row['Pathway_original_in_refV'].split('|'))-1

plotLengthDistributionAll()
plotLengthDistributionAllCrop()