# mamba create -n regen2 python=3.12.0
# mamba install pysam
# pip install pybedtools

import os, sys, random, pybedtools, statistics, argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import helper_functions

region_file_columns = ["chr", "start", "end", "region_id", "is_hit", "strand"]
result_file_columns = ['annotation', 'observed', 'fc', 'pval',
                       'hit_overlap', 'hit_fc', 'hit_pval',
                       'nonhit_overlap', 'nonhit_fc', 'nonhit_pval', 
                       'hit_overlap_per_exon', 'nonhit_overlap_per_exon', 'pval_per_exon']
CONFIG = {}

def update_config(args):
    CONFIG['REGION_FILE'] = args.region
    CONFIG['ANNOTATION_FILE'] = args.annotation
    CONFIG['OUT_DIR'] = args.outdir

    out_dir = args.outdir
    
    work_dir = os.path.join(out_dir, 'tmp')
    os.makedirs(work_dir, exist_ok = True)
    CONFIG['WORK_DIR'] = work_dir

    plot_dir = os.path.join(out_dir, 'plots')
    os.makedirs(plot_dir, exist_ok = True)
    CONFIG['PLOT_DIR'] = plot_dir

    annotation_file_base_name = os.path.basename(args.annotation)
    annotation_file_name = os.path.splitext(annotation_file_base_name)[0]
    out_dir = os.path.join(out_dir, annotation_file_name)
    os.makedirs(out_dir, exist_ok = True)
    CONFIG['ANNOTATION_DIR'] = out_dir

    CONFIG['LENGTH_TOLERENCE'] = args.length_tolerence
    CONFIG['SIMULATION_COUNT'] = args.simulation_count
    CONFIG['REGION_LENGTH_THRSHOLD'] = args.region_length_threshold

def intersect_region_annotation():
    region_file = CONFIG['REGION_FILE']
    annotation_file = CONFIG['ANNOTATION_FILE']
    
    region_bed = pybedtools.BedTool(region_file)
    annotation_bed = pybedtools.BedTool(annotation_file)

    # stranded intersection??
    intersection = region_bed.intersect(annotation_bed, wao=True)
    # intersection = region_exon_intersection.sort().merge(c='5,4,6', )

    annotation_file_base_name = os.path.basename(annotation_file)
    annotation_file_name = os.path.splitext(annotation_file_base_name)[0]
    intersection_file_name = 'regions_%s_intersection.tsv'%(annotation_file_name)
    CONFIG['INTERSECTION_FILE'] = os.path.join(CONFIG['ANNOTATION_DIR'], intersection_file_name)
    
    intersection.saveas(CONFIG['INTERSECTION_FILE'])

    intersectionDf = pd.read_csv(CONFIG['INTERSECTION_FILE'], sep='\t', header=None)
    
    # aggregate as counts if only intersted in number of peaks instead of intersection length
    # intersectionDf = intersectionDf.groupby(list(range(6))).agg({intersectionDf.columns[-1]:'count'}).reset_index()
    
    intersectionDf = intersectionDf.groupby(list(range(6))).agg({intersectionDf.columns[-1]:'sum'}).reset_index()

    intersection_file_name = 'regions_%s_intersection_grouped.csv'%(annotation_file_name)
    CONFIG['INTERSECTION_FILE_GROUPED'] = os.path.join(CONFIG['ANNOTATION_DIR'], intersection_file_name)

    intersectionDf.to_csv(CONFIG['INTERSECTION_FILE_GROUPED'], index=None)

# convert the input region_file into required format
def process_files(region_file, exon_file, annotation_file, out_dir):
    
    print(region_file, exon_file, annotation_file, out_dir)
    
    # process the region file?
    
    # intersect regions with exons
    
    region_bed = os.path.join(work_dir, 'regions.bed')
    print(region_bed)
    
    regions_df = pd.read_csv(region_file, sep='\t')
    regions_df.columns = region_file_columns
    # regions_df = regions_df[region_file_columns]
    regions_df.to_csv(region_bed, sep='\t', header=None, index=None)
    
    if len(exon_file)>0:
        region_bed = pybedtools.BedTool(region_bed)
        exon_bed = pybedtools.BedTool(exon_file)
        
        region_exon_intersection_bed = os.path.join(work_dir, 'regions_exon_intersect.bed')
        
        region_exon_intersection = region_bed.intersect(exon_bed)
        region_exon_intersection = region_exon_intersection.sort().merge(c='5,4,6', )
        region_exon_intersection.saveas(region_exon_intersection_bed)
        
    # process annotation file? seperate script
    # store the annotation file in cwd
    
def getRegionHash(regionList):
    regionList.sort()
    return '-'.join(regionList)

def getHitHash(hitList):
    hitList = [str(x) for x in hitList]
    return ''.join(hitList)

# returns a csv file with region id as first column and hit status as remaining columns
def generate_length_simulations(LengthTolerance):

    region_file = CONFIG['REGION_FILE']
    
    simulationFileName = 'length_simulation_%0.2f'%(LengthTolerance)
    simulationFile = os.path.join(CONFIG['WORK_DIR'], '%s.csv'%(simulationFileName))
    
    if os.path.exists(simulationFile):
        return simulationFile

    print('\nGenerating %s'%(simulationFile))
    
    regionsDf = pd.read_csv(region_file, sep='\t', header=None)
    regionsDf.columns = region_file_columns
    

    regionsDf['length'] = regionsDf['end'] - regionsDf['start']

    # length based filtering
    regionsDf = regionsDf[regionsDf['length']<=CONFIG['REGION_LENGTH_THRSHOLD']].reset_index(drop=True)
    
    hitRegionsDf = regionsDf[regionsDf['is_hit']==True]
    nonHitRegionsDf = regionsDf[regionsDf['is_hit']==False]
    
    regionListOriginal = list(regionsDf['region_id'])

    TotalHitLength = sum(hitRegionsDf['length'])
    LengthDiffThreshold = LengthTolerance * TotalHitLength

    print("%s * %s = %s"%(TotalHitLength, LengthTolerance, LengthDiffThreshold))

    hitLengthDistribution = os.path.join(CONFIG['PLOT_DIR'], '.'.join(['observed_length_distribution', 'png']))
    g = sns.kdeplot(hitRegionsDf, x="length", color='red')
    g = sns.kdeplot(nonHitRegionsDf, x="length")
    plt.savefig(hitLengthDistribution, bbox_inches='tight')
    plt.clf()

    # store the hash of each simulation so that they are not repeated 
    simulationHash = set()
    # simulationHash.add(getRegionHash(regionsDf[regionsDf['is_hit']]['region_id']))
    simulationHash.add(getHitHash(regionsDf['is_hit']))
    # random.seed(101)

    rowCount = len(regionsDf.index)
    rows = range(len(regionsDf['region_id']))
    cols = ['region_id']
    cols.extend([str(x) for x in range(CONFIG['SIMULATION_COUNT'])])
    simulationDf = pd.DataFrame(index=rows, columns=cols)
    simulationDf['region_id'] = regionsDf['region_id']
    
    i = 0
    unusedSimulations = 0
    # unusedSimulationsthreshold = 0.1*CONFIG['SIMULATION_COUNT']
    unusedSimulationsthreshold = 0.5*CONFIG['SIMULATION_COUNT']

    fig, axs = plt.subplots(ncols=2,nrows=1,figsize=(20,10),dpi=200)
    axs[0].set_title('Positive class')
    axs[1].set_title('Negative class')
    
    while i < CONFIG['SIMULATION_COUNT']:
        # print(i)
        currLen = 0
        saveFlag = True

        hitSet = set()
        regionList = []

        # pick a region as hit till the total hit length is within a range
        while True:
            r = random.randrange(0, rowCount)
            while r in hitSet:
                r = random.randrange(0, rowCount)

            hitSet.add(r)
            row = regionsDf.iloc[r]
            lr = row['length']
            currLen += lr
            regionId = row['region_id']
            regionList.append(regionId)

            lengthDiff = abs(currLen - TotalHitLength)
            if lengthDiff < LengthDiffThreshold:
                break
            elif currLen > TotalHitLength:
                print("\tSimulation %s dropped: total length (%s) exceeds the threshold"%(str(i), str(currLen)))
                saveFlag = False
                break
                
        # regionHash = getRegionHash(regionList)
        # if regionHash in simulationHash or not saveFlag:
        if not saveFlag:
            # print("Simulation %s dropped: already exists"%(str(i)))
            # print(regionHash)
            unusedSimulations += 1
        else:
            colName = str(i)
            simulationDf[colName] = 0
            
            simulationDf.loc[simulationDf['region_id'].isin(regionList), colName] = 1

            regionHash = getHitHash(simulationDf[colName])

            if regionHash in simulationHash:
                unusedSimulations += 1
            else:
                i += 1
                simulationHash.add(regionHash)
                
                hitDf = regionsDf[regionsDf['region_id'].isin(regionList)].copy()
                g1 = sns.kdeplot(hitDf, x="length", ax=axs[0])
                
                nonHitDf = regionsDf[~regionsDf['region_id'].isin(regionList)].copy()
                g2 = sns.kdeplot(nonHitDf, x="length", ax=axs[1])
            
        if unusedSimulations > unusedSimulationsthreshold:
            print(unusedSimulations, unusedSimulationsthreshold)
            warningMessage = 'Simulation limit exceeded (%d)! Only %d simulation completed out of %d'%(unusedSimulations, i, CONFIG['SIMULATION_COUNT'])
            print(warningMessage)
            simulationDf = simulationDf[simulationDf.columns[:i+1]].copy()
            break

    lengthDistribution = os.path.join(CONFIG['PLOT_DIR'], '%s_simulated_length_distribution.png'%(simulationFileName))
    plt.savefig(lengthDistribution, bbox_inches='tight')
    plt.clf()
    
    simulationDf.to_csv(simulationFile, index=None)
    return simulationFile

# returns a csv file with region id as first column and hit status as remaining columns
def generate_count_simulations(LengthTolerance, same_length=True):

    region_file = CONFIG['REGION_FILE']

    simulationFileName = 'both' if same_length else 'count'
    simulationFileName = '%s_simulation_%0.2f'%(simulationFileName, LengthTolerance)
    simulationFile = os.path.join(CONFIG['WORK_DIR'], '%s.csv'%(simulationFileName))
    
    if os.path.exists(simulationFile):
        return simulationFile

    print('\nGenerating %s'%(simulationFile))
    
    regionsDf = pd.read_csv(region_file, sep='\t', header=None)
    regionsDf.columns = region_file_columns

    regionsDf['length'] = regionsDf['end'] - regionsDf['start']

    # length based filtering
    regionsDf = regionsDf[regionsDf['length']<=CONFIG['REGION_LENGTH_THRSHOLD']].reset_index(drop=True)

    rowCount = len(regionsDf.index)
    
    hitRegionsDf = regionsDf[regionsDf['is_hit']==True]
    nonHitRegionsDf = regionsDf[regionsDf['is_hit']==False]
    
    isHitArrayOriginal = np.array(regionsDf['is_hit'])
    lengthArrayOriginal = np.array(regionsDf['length'])

    TotalHitLength = sum(hitRegionsDf['length'])
    LengthDiffThreshold = LengthTolerance * TotalHitLength

    # store the hash of each simulation so that they are not repeated 
    simulationHash = set()
    simulationHash.add(getHitHash(isHitArrayOriginal))
    # random.seed(101)

    rows = range(len(regionsDf['region_id']))
    cols = ['region_id']
    cols.extend([str(x) for x in range(CONFIG['SIMULATION_COUNT'])])
    simulationDf = pd.DataFrame(index=rows, columns=cols)
    simulationDf['region_id'] = regionsDf['region_id']
    
    i = 0
    unusedSimulations = 0
    # higher threshold is required for stringent simulation i.e. if both count and length need to be matched
    unusedSimulationsthreshold = 0.8*CONFIG['SIMULATION_COUNT'] if same_length else 0.1*CONFIG['SIMULATION_COUNT'] 

    fig, axs = plt.subplots(ncols=2,nrows=1,figsize=(20,10),dpi=200)
    axs[0].set_title('Positive class')
    axs[1].set_title('Negative class')
    
    while i < CONFIG['SIMULATION_COUNT']:
        # shuffle in place
        isHitArray = np.copy(isHitArrayOriginal)
        np.random.shuffle(isHitArray)
        
        currLen = np.sum(np.multiply(isHitArray, lengthArrayOriginal))
        lengthDiff = abs(currLen - TotalHitLength)
        
        hitHash = getHitHash(isHitArray)

        if hitHash in simulationHash or (lengthDiff > LengthDiffThreshold and same_length):
            # print("Simulation %s dropped: already exists"%(str(i)))
            # print(regionHash)
            unusedSimulations += 1
        else:
            colName = str(i)
            simulationDf[colName] = isHitArray
            regionList = simulationDf[simulationDf[colName]==1]['region_id']
            i += 1
            simulationHash.add(hitHash)
            
            hitDf = regionsDf[regionsDf['region_id'].isin(regionList)].copy()
            g1 = sns.kdeplot(hitDf, x="length", ax=axs[0])
            
            nonHitDf = regionsDf[~regionsDf['region_id'].isin(regionList)].copy()
            g2 = sns.kdeplot(nonHitDf, x="length", ax=axs[1])
            
        if unusedSimulations > unusedSimulationsthreshold:
            warningMessage = 'Simulation limit exceeded (%d)! Only %d simulation completed out of %d'%(unusedSimulations, i, CONFIG['SIMULATION_COUNT'])
            print(warningMessage)
            simulationDf = simulationDf[simulationDf.columns[:i+1]].copy()
            break

    lengthDistribution = os.path.join(CONFIG['PLOT_DIR'], '%s_simulated_length_distribution.png'%(simulationFileName))
    plt.savefig(lengthDistribution, bbox_inches='tight')
    plt.clf()
    
    simulationDf.to_csv(simulationFile, index=None)
    return simulationFile

def getMetrics(expectedList, observed, plotPath):

    print('Calculating empitical p value based on %d simulations'%len(expectedList))

    helper_functions.sns_set_my_style('white')
    expectedListMedian = statistics.median(expectedList)
    
    # p value is the probability of finding more extreme result than the observed ratio
    if observed>=expectedListMedian:
        p = sum([int(x>=observed) for x in expectedList])/len(expectedList)
    else:
        p = sum([int(x<=observed) for x in expectedList])/len(expectedList)

    g = sns.kdeplot(expectedList)
    xList = [expectedListMedian, observed]
    cList = ['blue', 'red']

    try:
        maxDensity = np.max(g.lines[0].get_ydata())
    except:
        maxDensity = 0.1
    
    plt.vlines(xList, 0, maxDensity, colors=cList, linestyles='dashdot')
    plt.xlabel('Hit/Non-hit overlap ratio')

    # annotate with fc and pvalue
    plt.text(1.05*expectedListMedian, 3*maxDensity/4, 'Expected = %.3f'%(expectedListMedian))
    plt.text(1.05*observed, maxDensity/2, 'Observed = %.3f'%(observed))
    plt.text(1.05*observed, maxDensity/4, 'p = %.3f'%(p))
    
    plt.savefig(plotPath, bbox_inches='tight')
    plt.clf()

    ratio = round(observed/expectedListMedian, 2) if expectedListMedian>0 else 0
    # ratio = round(observed/(expectedListMedian+1), 2)
    
    return ratio, p

def calculate_significance(simulationFile):
    intersection_file = CONFIG['INTERSECTION_FILE_GROUPED']
    elementName = os.path.basename(CONFIG['ANNOTATION_DIR'])
    intersectionDf = pd.read_csv(intersection_file)
    intersectionDf.columns = [int(x) for x in intersectionDf.columns]
    # print(intersectionDf.head())
    simulationFileDf = pd.read_csv(simulationFile)
    colList = simulationFileDf.columns[1:]
    
    simulationFileDf = pd.merge(simulationFileDf, intersectionDf, left_on='region_id', right_on=3)
    # print(simulationFileDf.head())

    simulationFileDf['l'] = simulationFileDf[2] - simulationFileDf[1]
    simulationFileDf['percent_overlap'] = 100*simulationFileDf[12]/simulationFileDf['l']
    
    # length based filtering
    # simulationFileDf = simulationFileDf[simulationFileDf['l']<=CONFIG['REGION_LENGTH_THRSHOLD']].copy()

    expectedRatioList = []
    expectedHitList = []
    expectedNonHitList = []

    subsetColList = [12, 'l']
    
    for colName in colList:
        
        hitDf = simulationFileDf[simulationFileDf[colName]==1][subsetColList]
        
        nonhitDf = simulationFileDf[simulationFileDf[colName]==0][subsetColList]

        # need to add psuedo count to avoid division by zero when calculating ratio
        
        hitOverlap = 100*(sum(hitDf[12]) + 1)/(sum(hitDf['l']) + 1)
        expectedHitList.append(round(hitOverlap, 2))
        
        nonhitOverlap = 100*(sum(nonhitDf[12]) + 1)/(sum(nonhitDf['l']) + 1)
        expectedNonHitList.append(round(nonhitOverlap, 2))
        
        overlapRatio = hitOverlap/nonhitOverlap
        expectedRatioList.append(round(overlapRatio, 2))

    file_base_name = os.path.basename(simulationFile)
    file_name = os.path.splitext(file_base_name)[0]

    # calculate observed ratio
    colName = 4
    
    hitDf = simulationFileDf[simulationFileDf[colName]==1]
        
    nonhitDf = simulationFileDf[simulationFileDf[colName]==0]

    hitOverlap = 100*(sum(hitDf[12]) + 1)/(sum(hitDf['l']) + 1)
    hitOverlap = round(hitOverlap, 2)
    ratioDistribution = os.path.join(CONFIG['ANNOTATION_DIR'], '%s_hit_overlap_distribution.png'%(file_name))
    hitRatio, hitPval = getMetrics(expectedHitList, hitOverlap, ratioDistribution)
    
    nonhitOverlap = 100*(sum(nonhitDf[12]) + 1)/(sum(nonhitDf['l']) + 1)
    nonhitOverlap = round(nonhitOverlap, 2)
    ratioDistribution = os.path.join(CONFIG['ANNOTATION_DIR'], '%s_nonhit_overlap_distribution.png'%(file_name))
    nonhitRatio, nonhitPval = getMetrics(expectedNonHitList, nonhitOverlap, ratioDistribution)
    
    overlapRatio = hitOverlap/nonhitOverlap if nonhitOverlap>0 else 0
    # overlapRatio = hitOverlap/(nonhitOverlap+1)
    overlapRatio = round(overlapRatio, 2)
    ratioDistribution = os.path.join(CONFIG['ANNOTATION_DIR'], '%s_overlap_ratio_distribution.png'%(file_name))
    rfc, rPval = getMetrics(expectedRatioList, overlapRatio, ratioDistribution)
    
    # helper_functions.custom_box_plot(simulationFileDf, colName, 'percent_overlap', ylabel='%s percent overlap'%(file_name), xlabel='Is exon hit?', yLimMin=-5, showViolin=False)
    pval = helper_functions.plot_hit_nonhit_cdf(simulationFileDf, colName, 'percent_overlap', titleStr=elementName)
    hitExonMedianOverlap = np.median(hitDf['percent_overlap'])
    nonhitExonMedianOverlap = np.median(nonhitDf['percent_overlap'])
    
    plotPath = os.path.join(CONFIG['ANNOTATION_DIR'], '%s_overlap_distribution.png'%(file_name))
    plt.savefig(plotPath, bbox_inches='tight')
    plt.clf()

    # save the results to a file
    writeList = [os.path.splitext(os.path.basename(CONFIG['ANNOTATION_FILE']))[0],
                 overlapRatio, rfc, rPval,
                 hitOverlap, hitRatio, hitPval,
                 nonhitOverlap, nonhitRatio, nonhitPval, hitExonMedianOverlap, nonhitExonMedianOverlap, pval]
    
    writeList = [str(x) for x in writeList]

    regionName = os.path.splitext(os.path.basename(CONFIG['REGION_FILE']))[0]
    simulationName = os.path.splitext(os.path.basename(simulationFile))[0]
    savePath = os.path.join(CONFIG['OUT_DIR'], '%s_enrichment_%s.tsv'%(regionName, simulationName))

    if os.path.exists(savePath):
        with open(savePath, 'a') as f:
            f.write('\t'.join(writeList) + '\n')
    else:
        with open(savePath, 'w') as f:
            f.write('\t'.join(result_file_columns) + '\n')
            f.write('\t'.join(writeList) + '\n')

    # tmpDf = pd.read_csv(savePath, sep='\t')
    # tmpDf = tmpDf.sort_values(by=['fc'], ascending=False)
    # tmpDf.to_csv(savePath, sep='\t', index=None)
    
    print('Results saved to %s'%(savePath)) 

def parse_args():
    parser = argparse.ArgumentParser(
        description="REGEN (REgion-based Genomic ENrichment): tool for evaluating enrichment of annotations within genomic regions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required argument
    parser.add_argument(
        "--region",
        required=True,
        help="Path to region file, should be in bed6 format with 4th column as region_id and 5th column as is_hit",
    )

    # Required argument
    parser.add_argument(
        "--annotation",
        required=True,
        help="Path to annotation file, should be sorted and merged",
    )

    # Required argument
    parser.add_argument(
        "--outdir",
        required=True,
        help="Path of output folder",
    )

    parser.add_argument(
        "--mode",
        type=str,
        default='LENGTH',
        help="maintain either LENGTH, COUNT or BOTH during simulations",
    )

    parser.add_argument(
        "--length_tolerence",
        type=float,
        default=0.05,
        help="Threshold of length for the simulated regions as fraction of the actual length of the regions",
    )

    parser.add_argument(
        "--region_length_threshold",
        type=int,
        default=2000,
        help="Only consider regions which are shorter than this threshold",
    )

    parser.add_argument(
        "--simulation_count",
        type=int,
        default=10000,
        help="Count of simulated regions to calcululate sampling distribution of the overlap",
    )
    
    return parser.parse_args()

if __name__=='__main__':

    args = parse_args()
    update_config(args)
    
    intersect_region_annotation()

    print(CONFIG)

    print('\nCurrent legth tolerance: %f'%(args.length_tolerence))

    if args.mode=='LENGTH':
        simulationFile = generate_length_simulations(args.length_tolerence)
    else:
        simulationFile = generate_length_simulations(args.length_tolerence, same_length=args.mode=='BOTH')
    print(simulationFile)
    
    calculate_significance(simulationFile)