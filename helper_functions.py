from ast import literal_eval
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager
import seaborn as sns

from scipy.stats import mannwhitneyu, ttest_ind
from itertools import combinations

# GRAPHICAL OPTIONS
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16
def sns_set_my_style(s, font_family=''):#white,whitegrid
    # sns.set_style(s)
    sns.set_theme(style=s, palette='colorblind')
    if 'grid' in s:
        plt.rc('grid', alpha=0.5)
        plt.rc('grid', linestyle='--')

    # matplotlib.font_manager.fontManager.addfont('/home/people/22211214/scratch/fonts/arial/ARIAL.TTF')
    font_family = font_family if font_family else 'Arial'
    plt.rc('font', size=SMALL_SIZE, family=font_family)
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    plt.rc('figure', dpi=300)
    plt.rc('savefig', bbox='tight')

# coolwarm
def custom_box_plot(plotData, xCol, yCol, savePath='', dropZeroes=False, xlabel='', ylabel='', ax='', title='', showViolin=True, showPoints=False, hideNonSignificantPvalue=False, showPvalues=True, yLimMax='', yLimMin='', palette='colorblind', figSize='', roundDigit=4, showGroupCount=True):
    
    if not ax:
        fig,ax = plt.subplots(ncols=1,nrows=1,figsize=figSize) if figSize else plt.subplots(ncols=1,nrows=1)

    plotData = plotData.dropna(subset=[xCol, yCol])
    
    if dropZeroes:
        plotData = plotData[plotData[yCol]!=0]
    
#     if showViolin:
#         # violin_plot = sns.violinplot(x=xCol, y=yCol, data=plotData, dodge=False, scale="width", inner=None, cut=0.5, ax=ax)
#         violin_plot = sns.violinplot(x=xCol, y=yCol, data=plotData, dodge=False, density_norm="width", inner=None, cut=0.5, ax=ax)
    
#     box_plot = sns.boxplot(x=xCol, y=yCol, data=plotData, saturation=1, showfliers=True,
#             width=0.3, boxprops={'zorder': 3, 'facecolor': 'none'}, ax=ax)
    
    if showViolin:
        # violin_plot = sns.violinplot(x=xCol, y=yCol, data=plotData, dodge=False, scale="width", inner=None, cut=0.5, ax=ax)
        violin_plot = sns.violinplot(x=xCol, y=yCol, data=plotData, dodge=False, density_norm="width", inner=None, cut=0.5, ax=ax, hue=xCol, palette=palette)
        box_plot = sns.boxplot(x=xCol, y=yCol, data=plotData, saturation=1, showfliers=False, width=0.3, ax=ax, boxprops={'zorder': 3, 'facecolor': 'none'}, hue=xCol, palette=palette)
    elif showPoints:
        box_plot = sns.boxplot(x=xCol, y=yCol, data=plotData, saturation=1, showfliers=False, width=0.3, ax=ax, medianprops={"linewidth": 3}, hue=xCol, palette=palette)
        sns.stripplot(data=plotData, x=xCol, y=yCol, dodge=True, ax=ax, alpha=0.5, color='black')
    else:
        box_plot = sns.boxplot(x=xCol, y=yCol, data=plotData, saturation=1, showfliers=True, width=0.3, ax=ax, medianprops={"linewidth": 3}, hue=xCol, palette=palette)
    
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    
    medians = plotData.groupby([xCol])[yCol].median()
    # rowCounts = psiDfFiltered.groupby(['isHit'])['psi'].shape[0]
    vertical_offset = plotData[yCol].median() * 0.1 # offset from median for display
    # vertical_offset = 0.05
    
    # xticks = box_plot.get_xticks()
    xticks = medians.index
    xcount = len(xticks)
    
    # print(medians)
    # print(xticks)

    yMax = plotData[yCol].max()
    yMin = plotData[yCol].min()
    
    yLimMax = yLimMax if yLimMax else yMax
    yLimMin = yLimMin if yLimMin else yMin
    
    # plt.ylim(-0.05*yMax, 1.4*yMax)
    ax.set_ylim([1.1 * yLimMin, 1.6 * yLimMax])
    # ax.set_ylim([-0.05*yMax, (1.45 + 0.1*xcount) * yMax])
    
    for xtick in medians.index:
        nonZeroCount = plotData[(plotData[xCol]==xtick) & (plotData[yCol]!=0)].shape[0]
        totalCount = plotData[plotData[xCol]==xtick].shape[0]
        nonZeroPerc = round(100*nonZeroCount/totalCount, 2) if totalCount>0 else 0
        if xtick in medians:
            # textStr =  str(round(medians[xtick], 4)) + '\n' + 'n = ' + str(totalCount) + '\n' + 'Non-zero % = ' + str(nonZeroPerc) 
            textStr =  str(round(medians[xtick], roundDigit)) + '\n' + 'n = ' + str(totalCount) if showGroupCount else str(round(medians[xtick], roundDigit))
            # medians[xtick] + vertical_offset
            # box_plot.text(xtick, 1.1*yMax, textStr, horizontalalignment='center', size='x-small',color='b',weight='semibold')
            box_plot.text(xtick, 1.15*yMax, textStr, horizontalalignment='center',color='black',weight='semibold')
    
    if showPvalues:
        pairs = list(combinations(xticks, 2))  # Generate all pairwise combinations
        p_values = {}

        for pair in pairs:
            group1 = plotData[plotData[xCol]==pair[0]][yCol]
            group2 = plotData[plotData[xCol]==pair[1]][yCol]

            # group1 = np.array(group1) / np.max(group1)
            # group2 = np.array(group2) / np.max(group2)

            if len(group1)>0 and len(group2)>0:
                stat, p_val = mannwhitneyu(group1, group2)
                # Welch's t-test
                # stat, p_val = ttest_ind(group1, group2, equal_var=False)
                p_values[pair] = p_val
                # print(pair, p_val)

        step = (yLimMax * 0.2)  # Vertical spacing between annotations

        for i, ((cat1, cat2), p_val) in enumerate(p_values.items()):
            # print(p_val)
            pvalTextColor = 'black'
            if hideNonSignificantPvalue and p_val > 0.05:
                continue
            if hideNonSignificantPvalue:
                pvalTextColor = 'black'
            else:
                pvalTextColor = 'black' if p_val > 0.05 else 'red'
            x1, x2 = xticks.tolist().index(cat1), xticks.tolist().index(cat2)
            # y, h = yMax + (i + 3) * step, step * 0.5
            y, h = yMax + (i + 3) * step, step * 0.3
            ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, color="black")
            box_plot.text((x1 + x2) * 0.5, y + h, f"p = {p_val:.3f}", ha="center", va="bottom", color=pvalTextColor)
    
    axlegend = ax.get_legend()
    if axlegend:
        axlegend.remove()
    
    plt.tight_layout()
    
    if savePath:
        plt.savefig(savePath, bbox_inches='tight')

def plot_hit_nonhit_cdf(data, hitColName, colName, ax='', xStr='percent overlap', savePath='', titleStr=''):

    if not ax:
        fig,ax = plt.subplots(ncols=1,nrows=1)
        
    cmap = {0:'r', 1:'g'}
    dataList = []
    for h in [0, 1]:
        pdata = data[data[hitColName]==h][colName]
        dataList.append(pdata)
        # pdata = np.log2(pdata+1)
        # cdfPlot = sns.ecdfplot(pdata, ax=ax, color=cmap[h], label='is hit = %d'%(h))

    sns.ecdfplot(data, x=colName, ax=ax, palette=cmap, hue=hitColName)
    stat, p_val = mannwhitneyu(dataList[0], dataList[1])
    ax.text(0.3, 0.85, 'p=%.3f'%(p_val), transform=ax.transAxes)
    ax.set_ylabel('Cumulative Probability')
    ax.set_xlabel(xStr)
    ax.set_title(titleStr)

    plt.tight_layout()
    
    if savePath:
        plt.savefig(savePath, bbox_inches='tight')

    return p_val

# Plots the PDF and CDF of a given dataset on the same plot, with vertical lines for quartiles
def plot_distribution(data, log2Transform=False, pseudoCount='', title="PDF and CDF with Quartiles", xLabel='', ax='', savePath=''):
    
    if log2Transform:
        data = np.log2(data+pseudoCount) if pseudoCount else np.log2(data)
    
    # Create a figure and a primary axis
    if not ax:
        # fig, ax1 = plt.subplots()
        fig,ax1 = plt.subplots(ncols=1,nrows=1)
    else:
        ax1 = ax

    # Plot the PDF on the primary axis (left)
    # ax1.hist(data, bins=30, density=True, alpha=0.6, color='b', label='PDF')
    # sns.displot(exonDf[exonDf['l']<5000], x="l", kind='kde')
    pdfPlot = sns.kdeplot(data, ax=ax1)
    ax1.set_xlabel(xLabel)
    ax1.set_ylabel('Probability Density', color='b')
    ax1.tick_params(axis='y', labelcolor='b')

    # Create a secondary axis for the CDF (right)
    ax2 = ax1.twinx()

    # Plot the CDF on the secondary axis
    # sorted_data = np.sort(data)
    # cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
    # ax2.plot(sorted_data, cdf, color='r', label='CDF')
    cdfPlot = sns.ecdfplot(data, ax=ax2, color='r')
    ax2.set_ylabel('Cumulative Probability', color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    
    xcor = np.min(data)
    xcorLabel = np.exp2(xcor) if log2Transform else xcor
    label = 'Min = \n' + str(round(xcorLabel, 2))
    ax2.text(xcor, 0.1, label, horizontalalignment='right',color='black',weight='semibold')
    
    xcor = np.max(data)
    xcorLabel = np.exp2(xcor) if log2Transform else xcor
    label = 'Max = \n' + str(round(xcorLabel, 2))
    ax2.text(xcor, 0.1, label, horizontalalignment='right',color='black',weight='semibold')
    
    xcor = np.max(data)
    label = 'Count = ' + str(len(data))
    ax2.text(xcor, 0.6, label, horizontalalignment='right',color='black',weight='semibold')
    
    # Add title and legend
    ax2.set_title(title)
    # fig.legend(loc="upper left", bbox_to_anchor=(0.15, 0.9))
    
    # Calculate quartiles
    quartiles = np.percentile(data, [25, 50, 75])
    quartile_labels = ['25th Percentile (Q1)', '50th Percentile (Median)', '75th Percentile (Q3)']
    quartile_labels = ['Q1', 'Q2', 'Q3']
    
    # Draw vertical lines for each quartile
    i = 0.05
    for q, label in zip(quartiles, quartile_labels):
        qLabel = np.exp2(q) if log2Transform else q
        textStr = '%s = %s'%(label, str(round(qLabel, 2)))
        ax1.axvline(q, color='k', linestyle='--', alpha=0.7, label=textStr)
        ax2.text(q, 0.9, label, horizontalalignment='right',color='black',weight='semibold', rotation=20)
        ax2.text(xcor, 0.6-i, textStr, horizontalalignment='right',color='black',weight='semibold')
        i += 0.05

    plt.tight_layout()
    
    if savePath:
        plt.savefig(savePath)