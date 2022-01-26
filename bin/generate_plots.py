#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import math

argparser = argparse.ArgumentParser(description = 'Generates plots from the summary file.')
argparser.add_argument('-s', '--summary', metavar = 'file', dest = 'summary_filename', required = True, help = 'Name of the summary file.')
argparser.add_argument('-pca1', '--pca-1000g', metavar = 'file', dest = 'pca_1000g_filename', required = True, help = 'File with the PCs from 1000G.')
argparser.add_argument('-pop1', '--populations-1000g', metavar = 'file', dest = 'pop_1000g_filename', required = True, help = 'File with the populations from 1000G.')
argparser.add_argument('-pca2', '--pca-hgdp', metavar = 'file', dest = 'pca_hgdp_filename', required = False, help = 'File with the PCs from HGDP.')
argparser.add_argument('-pop2', '--populations-hgdp', metavar = 'file', dest = 'pop_hgdp_filename', required = False, help = 'File with the populations from HGDP.')
argparser.add_argument('-sex', '--reported-sex', metavar = 'file', dest = 'reported_sex_filename', required = True, help = 'Reported sex file. Tab-delimited. Header: NAME, SEX. Coding: F - female, M - male, NA - missing.')


colors_1000g = {
    'EUR': '#018ead',
    'EAS': '#778500',
    'SAS': '#c44cfd',
    'AMR': '#710027',
    'AFR': '#ffd845'
}


colors_hgdp = {
    'Africa': '#f2c2ff',
    'Middle_Est':  '#4e2100',
    'Central_South_Asia':  '#d8ffdf',
    'Est_Asia': '#b35200',
    'Oceania': '#001991',
    'Europe': '#284600',
    'America': '#ff9081'
}


def plot_average_depth(summary, output_filename):
    fig = plt.figure(figsize=(5, 5), dpi = 300)
    ax = fig.add_subplot(1, 1, 1)
    ax.violinplot(summary.AUTO_AVG_DP, showmedians = False, showmeans = False, showextrema = False)
    ax.boxplot(summary.AUTO_AVG_DP)
    ax.set_ylim([0, None])
    ax.set_ylabel('Average depth in autosomal chromosomes')
    ax.set_xticks([])
    ax.tick_params(axis = 'both', direction = 'in')
    ax.grid(linestyle = '--', linewidth = 0.5)
    plt.savefig(output_filename, bbox_inches='tight')


def plot_coverage(summary, output_filename):
    fig = plt.figure(figsize=(8, 5), dpi = 300)
    dps = [1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]
    data = []
    for dp in dps:
        data.append(summary[f'AUTO_PRCNT_DP{dp}'].values)
    ax = fig.add_subplot(1, 1, 1)
    ax.violinplot(data, showmedians = False, showmeans = False, showextrema = False)
    ax.boxplot(data, widths = 0.3, flierprops = {'markersize': 3})
    ax.tick_params(axis = 'both', direction = 'in')
    ax.grid(linestyle = '--', linewidth = 0.5)
    ax.set_ylim([0, 100])
    ax.set_yticks(range(0, 110, 10))
    ax.set_xticklabels(dps)
    ax.set_xlabel('Minimal read depth')
    ax.set_ylabel('Coverage (%)')
    plt.tight_layout()
    plt.savefig(output_filename, bbox_inches='tight')


def plot_contamination(summary, output_filename):
    fig = plt.figure(figsize=(4, 5), dpi = 300)
    ax = fig.add_subplot(1, 1, 1)
    data = [summary['1000G_FREEMIX'] * 100]
    if 'HGDP_FREEMIX' in summary.columns:
        data.append(summary['HGDP_FREEMIX'] * 100)
        ax.violinplot(data, showmedians = False, showmeans = False, showextrema = False)
    else:
        data.append([])
    ax.boxplot(data, widths = 0.3, flierprops = {'markersize': 3})
    ax.tick_params(axis = 'both', direction = 'in')
    ax.grid(linestyle = '--', linewidth = 0.5)
    ax.set_xticklabels(['1000G', 'HGDP'])
    ax.set_xlabel('Reference panel')
    ax.set_ylabel('Contamination (%)')
    plt.tight_layout()
    plt.savefig(output_filename, bbox_inches='tight')


def tabulate_contamination(summary, output_filename):
    contamination_table = []
    for low, high in [(-1, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.05), (0.05, 0.07), (0.07, 0.10), (0.10, 1.0)]:
        entry = {'% Contamination': f'{low * 100 if low > 0 else 0:.0f}-{high * 100:.0f}%'}
        entry['N samples (using 1000G)'] = len(summary[(summary['1000G_FREEMIX'] > low) & (summary['1000G_FREEMIX'] <= high)])
        if 'HGDP_FREEMIX' in summary.columns:
            entry['N samples (using HGDP)'] = len(summary[(summary['HGDP_FREEMIX'] > low) & (summary['HGDP_FREEMIX'] <= high)])
        contamination_table.append(entry)
    contamination_table = pd.DataFrame(contamination_table)    

    fig = plt.figure(figsize=(7, 2), dpi = 300)
    ax = plt.subplot(111, frame_on = False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    if 'HGDP_FREEMIX' in summary.columns:
        ax.table(cellText = contamination_table.values, colLabels = contamination_table.columns, loc = "center", colColours = ['lightgrey'] * 3, fontsize = 14)
    else:
        ax.table(cellText = contamination_table.values, colLabels = contamination_table.columns, loc = "center", colColours = ['lightgrey'] * 2, fontsize = 14)
    plt.savefig(output_filename)


def plot_pca_projection_all(summary, reference, reference_colors, reference_name, output_filename):
    fig = plt.figure(figsize=(10, 10), dpi = 300)

    for i, (pc1, pc2) in enumerate([('PC1', 'PC2'), ('PC1', 'PC3'), ('PC2', 'PC3'), ('PC2', 'PC4')], 1):
        ax = fig.add_subplot(2, 2, i)
        for pop, color in reference_colors.items():
            ax.scatter(reference[reference.POP == pop][pc1], reference[reference.POP == pop][pc2], facecolors = color, edgecolors = color, alpha = 0.3, label = pop)
        ax.scatter(summary[f'{reference_name}_{pc1}_INTENDED'], summary[f'{reference_name}_{pc2}_INTENDED'], facecolor = 'none', edgecolors = 'black', alpha = 1, label = 'Study')
        ax.tick_params(axis = 'both', direction = 'in')
        ax.grid(linestyle = '--', linewidth = 0.5)
        ax.set_xlabel(pc1)
        ax.set_ylabel(pc2)
        ax.legend()
    
    fig.suptitle(f"PCA projection to {reference_name}")
    plt.tight_layout()
    plt.savefig(output_filename, bbox_inches='tight')


def plot_pca_projection_most_contaminated(summary, reference, reference_colors, reference_name, output_filename):
    fig = plt.figure(figsize=(10, 10), dpi = 300)

    highest_contamination = summary[summary[f'{reference_name}_FREEMIX'] == summary[f'{reference_name}_FREEMIX'].max()]

    for i, (pc1, pc2) in enumerate([('PC1', 'PC2'), ('PC1', 'PC3'), ('PC2', 'PC3'), ('PC2', 'PC4')], 1):
        ax = fig.add_subplot(2, 2, i)
        for pop, color in reference_colors.items():
            ax.scatter(reference[reference.POP == pop][pc1], reference[reference.POP == pop][pc2], facecolors = color, edgecolors = color, alpha = 0.3, label = pop)
        ax.scatter(highest_contamination[f'{reference_name}_{pc1}_INTENDED'], highest_contamination[f'{reference_name}_{pc2}_INTENDED'], facecolor = 'black', edgecolors = 'black', alpha = 1, label = 'Intended')
        ax.scatter(highest_contamination[f'{reference_name}_{pc1}_CONTAMINATING'], highest_contamination[f'{reference_name}_{pc2}_CONTAMINATING'], facecolor = 'red', edgecolors = 'red', alpha = 1, label = 'Contaminating')
        ax.tick_params(axis = 'both', direction = 'in')
        ax.grid(linestyle = '--', linewidth = 0.5)
        ax.set_xlabel(pc1)
        ax.set_ylabel(pc2)
        ax.legend()
    
    contamination_value = highest_contamination[f'{reference_name}_FREEMIX'].values[0]
    fig.suptitle(f"PCA projection to {reference_name}: sample with the highest contamination of {contamination_value * 100:.1f}%")
    plt.tight_layout()
    plt.savefig(output_filename, bbox_inches='tight')


def plot_reported_sex(summary, summary_with_sex, output_filename):
    fig = plt.figure(figsize=(10, 10), dpi = 300)

    ax1 = fig.add_subplot(2, 2, 1)
    ax1.scatter(summary.X_NORM_AVG_DP, summary.Y_NORM_AVG_DP, facecolors = 'none', edgecolors = 'black', alpha = 1)
    ax1.set_ylabel('Normalized average depth on chromosome Y')
    ax1.tick_params(axis = 'both', direction = 'in')
    ax1.grid(linestyle = '--', linewidth = 0.5)
    ax1.legend(labels = [f'All (N={len(summary)})'])

    ax2 = fig.add_subplot(2, 2, 2, sharey = ax1, sharex = ax1)
    ax2.scatter(summary_with_sex[summary_with_sex.SEX.isnull()].X_NORM_AVG_DP,
            summary_with_sex[summary_with_sex.SEX.isnull()].Y_NORM_AVG_DP, facecolors = 'none', edgecolors = 'tab:gray', alpha = 1)
    ax2.tick_params(axis = 'both', direction = 'in')
    ax2.grid(linestyle = '--', linewidth = 0.5)
    ax2.legend(title = 'Reported sex:', labels = [f'NA (N={len(summary_with_sex[summary_with_sex.SEX.isnull()])})'])

    ax3 = fig.add_subplot(2, 2, 3, sharey = ax1, sharex = ax1)
    ax3.scatter(summary_with_sex[(~summary_with_sex.SEX.isnull()) & (summary_with_sex.SEX.str.startswith('F'))].X_NORM_AVG_DP,
            summary_with_sex[(~summary_with_sex.SEX.isnull()) & (summary_with_sex.SEX.str.startswith('F'))].Y_NORM_AVG_DP, facecolors = 'none', edgecolors = 'tab:orange', alpha = 1)
    ax3.tick_params(axis = 'both', direction = 'in')
    ax3.grid(linestyle = '--', linewidth = 0.5)
    ax3.set_ylabel('Normalized average depth on chromosome Y')
    ax3.set_xlabel('Normalized average depth on chromosome X')
    ax3.legend(title = 'Reported sex:', labels = [f'Female (N={len(summary_with_sex[(~summary_with_sex.SEX.isnull()) & (summary_with_sex.SEX.str.startswith("F"))])})'])

    ax4 = fig.add_subplot(2, 2, 4, sharey = ax1, sharex = ax1)
    ax4.scatter(summary_with_sex[(~summary_with_sex.SEX.isnull()) & (summary_with_sex.SEX.str.startswith('M'))].X_NORM_AVG_DP,
            summary_with_sex[(~summary_with_sex.SEX.isnull()) & (summary_with_sex.SEX.str.startswith('M'))].Y_NORM_AVG_DP, facecolors = 'none', edgecolors = 'tab:blue', alpha = 1)
    ax4.tick_params(axis = 'both', direction = 'in')
    ax4.grid(linestyle = '--', linewidth = 0.5)
    ax4.set_xlabel('Normalized average depth on chromosome X')
    ax4.legend(title = 'Reported sex:', labels = [f'Male (N={len(summary_with_sex[(~summary_with_sex.SEX.isnull()) & (summary_with_sex.SEX.str.startswith("M"))])})'])

    fig.suptitle("Normalized read depth in X/Y chromosomes vs reported sex")
    plt.tight_layout()
    plt.savefig(output_filename, bbox_inches='tight')


def infer_genetic_ancestry(summary, reference, reference_name):
    counts = {
        6: Counter(),
        7: Counter(),
        8: Counter(),
        9: Counter(),
        10: Counter()
    }

    def euclidean_distance(v1, v2):
        return math.sqrt(sum((x - y)**2 for x, y in zip(v1, v2)))

    for index, study_sample in summary.iterrows():
        study_coords = [study_sample[f'{reference_name}_PC1_INTENDED'], study_sample[f'{reference_name}_PC2_INTENDED'], study_sample[f'{reference_name}_PC3_INTENDED'], study_sample[f'{reference_name}_PC4_INTENDED']]
        distances = reference[['PC1', 'PC2', 'PC3', 'PC4', 'POP']].copy()
        distances['DIST'] = distances.apply(lambda row: euclidean_distance(study_coords, [row['PC1'], row['PC2'], row['PC3'], row['PC4']]), axis = 1)
        pop, n = Counter(distances.sort_values(by = 'DIST')['POP'][:10]).most_common(1)[0]
        if n >= 6:
            for i in range(6, n + 1):
                counts[i][pop] += 1
            for i in range(n + 1, 11):
                counts[i]['Undefined'] += 1
        else:
            for i in range(6, 11):
                counts[i]['Undefined'] += 1

    return counts


def tabulate_inferred_genetic_ancestry(inferred_ancestries, reference_colors, reference_name, output_filename):
    inferred_ancestries_table = []
    for i in range(6, 11):
        entry = { 'k-NN': i }
        for pop in list(reference_colors.keys()) + ['Undefined']:
            entry[pop] = inferred_ancestries[i][pop]
        inferred_ancestries_table.append(entry)
    inferred_ancestries_table = pd.DataFrame(inferred_ancestries_table)

    fig = plt.figure(figsize=(13, 3), dpi = 300)
    ax = plt.subplot(111, frame_on = False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    table = ax.table(cellText = inferred_ancestries_table.values, colLabels = inferred_ancestries_table.columns, loc = "center", colColours = ['lightgrey'] * (len(reference_colors) + 2), fontsize = 14)
    ax.set_title(f'Inferred genetic ancestry based on k-nearest neighbours (k-NN) in {reference_name}', y = 0.75)
    table.auto_set_font_size(False)
    plt.tight_layout()
    plt.savefig(output_filename, bbox_inches='tight')


if __name__ == '__main__':
    args = argparser.parse_args()

    summary = pd.read_csv(args.summary_filename, header = 0, sep = ' ')

    plot_average_depth(summary, 'average_depth.jpeg')
    plot_coverage(summary, 'coverage.jpeg')
    plot_contamination(summary, 'contamination_boxplot.jpeg')
    tabulate_contamination(summary, 'contamination_table.jpeg')

    pca_1000g = pd.read_csv(args.pca_1000g_filename, sep = '\t', header = None, names = ['SAMPLE', 'PC1', 'PC2', 'PC3', 'PC4'], usecols = [0, 1, 2, 3, 4])
    pop_1000g = pd.read_csv(args.pop_1000g_filename, sep = '\t')
    assert len(pca_1000g) == len(pop_1000g), f'{len(pca_hgdp)} != {len(samples_hgdp)}'
    samples_1000g = pca_1000g.merge(pop_1000g)
    assert len(samples_1000g) == len(pca_1000g), f'{len(pca_hgdp)} != {len(samples_hgdp)}'
    plot_pca_projection_all(summary, samples_1000g, colors_1000g, '1000G', '1000G_pca.jpeg')
    plot_pca_projection_most_contaminated(summary, samples_1000g, colors_1000g, '1000G', '1000G_pca_most_contaminated.jpeg')

    inferred_ancestry = infer_genetic_ancestry(summary, samples_1000g, '1000G')
    tabulate_inferred_genetic_ancestry(inferred_ancestry, colors_1000g, '1000G', '1000G_inferred_ancestry_table.jpeg')
 
    if args.pop_hgdp_filename is not None and args.pca_hgdp_filename is not None:
        pca_hgdp = pd.read_csv(args.pca_hgdp_filename, sep = '\t', header = None, names = ['SAMPLE', 'PC1', 'PC2', 'PC3', 'PC4'], usecols = [0, 1, 2, 3, 4])
        pop_hgdp = pd.read_csv(args.pop_hgdp_filename, sep = '\t')
        samples_hgdp = pca_hgdp.merge(pop_hgdp)
        assert len(samples_hgdp) == len(pca_hgdp), f'{len(pca_hgdp)} != {len(samples_hgdp)}'
        plot_pca_projection_all(summary, samples_hgdp, colors_hgdp, 'HGDP', 'HGDP_pca.jpeg')
        plot_pca_projection_most_contaminated(summary, samples_hgdp, colors_hgdp, 'HGDP', 'HGDP_pca_most_contaminated.jpeg')
        inferred_ancestry = infer_genetic_ancestry(summary, samples_hgdp, 'HGDP')
        tabulate_inferred_genetic_ancestry(inferred_ancestry, colors_hgdp, 'HGDP', 'HGDP_inferred_ancestry_table.jpeg')

    reported_sex = pd.read_csv(args.reported_sex_filename, sep = '\t')
    summary_with_sex = summary.merge(reported_sex, on = 'NAME', how = 'left')
    plot_reported_sex(summary, summary_with_sex, 'XY_depth_vs_reported_sex.jpeg')

   

