#!/usr/bin/env python3

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import pandas as pd
import inspect 

plt.rcParams['svg.fonttype'] = 'none'

#=======================#
# Define distance files #
#=======================#

# Each of theses is the output of `bedtools closest` between the labeled
# novel regions and the named variant dataset

gwas_distance_file = 'closest_gwas_catalog_to_each_novel_region.tsv'
clinvar_distance_file = 'closest_clinvar_to_each_novel_region.tsv'
pathogenic_clinvar_distance_file = 'closest_pathogenic_clinvar_to_each_novel_region.tsv'


#==================#
# Define functions #
#==================#

def get_distances(file_path):
	out_distances = []
	out_labels = []
	for line in open(file_path):
		fields = line.rstrip('\n').split('\t')
		out_distances.append(int(fields[-1])/1000000) # Converts distances to Mbp
		out_labels.append(fields[3])
	out_df = pd.DataFrame({'distance': out_distances, 'label':out_labels})
	return out_df

def get_zoomed_distances(file_path, cutoff):
	out_distances = []
	out_labels = []
	for line in open(file_path):
		fields = line.rstrip('\n').split('\t')
		if int(fields[-1]) <= cutoff:
			out_distances.append(int(fields[-1])/1000000)
			out_labels.append(fields[3])
	out_df = pd.DataFrame({'distance': out_distances, 'label':out_labels})
	return out_df

def move_legend(ax, new_loc, **kws):
    old_legend = ax.legend_
    handles = old_legend.legendHandles
    labels = [t.get_text() for t in old_legend.get_texts()]
    title = old_legend.get_title().get_text()
    ax.legend(handles, labels, loc=new_loc, title=title, **kws)


#==================#
# Generate Figures #
#==================#

gwas_distances = get_distances(gwas_distance_file)
clinvar_distances = get_distances(clinvar_distance_file)
pathogenic_clinvar_distances = get_distances(pathogenic_clinvar_distance_file)

zoomed_cutoff = 1000000

gwas_zoomed_distances = get_zoomed_distances(gwas_distance_file, zoomed_cutoff)
clinvar_zoomed_distances = get_zoomed_distances(clinvar_distance_file, zoomed_cutoff)
pathogenic_clinvar_zoomed_distances = get_zoomed_distances(pathogenic_clinvar_distance_file, zoomed_cutoff)


sns.set_theme()
sns.set_style('white')

pal = ['#7995C4', '#592318', '#898A00']

has_legend = False

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True, figsize=(10,10))
sns.histplot(data=gwas_distances, ax=ax1, stat="count", bins=[x/10 for x in range(0,210,5)], multiple='stack', x='distance', kde=False, hue='label', element='bars', linewidth=0.5, palette=pal, legend=True)
sns.histplot(data=clinvar_distances, ax=ax2, stat="count", bins=[x/10 for x in range(0,210,5)], multiple='stack', x='distance', kde=False, hue='label', element='bars', linewidth=0.5, palette=pal, legend=has_legend)
sns.histplot(data=pathogenic_clinvar_distances, ax=ax3, stat="count", bins=[x/10 for x in range(0,210,5)], multiple='stack', x='distance', kde=False, hue='label', element='bars', palette=pal, linewidth=0.5, legend=has_legend)
for ax in [ax1, ax2, ax3]:
	ax.set_xlabel(None)
ax1.legend_.set_title(None)
move_legend(ax1, 'lower right')
ax1ins = inset_axes(ax1, width="60%", height="40%")
ax2ins = inset_axes(ax2, width="60%", height="40%")
ax3ins = inset_axes(ax3, width="60%", height="40%")
sns.histplot(data=gwas_zoomed_distances, ax=ax1ins, stat='count', bins=[x/50 for x in range(0,51)], multiple='stack', x='distance', kde=False, hue='label', element='bars', linewidth=0.5, palette=pal, legend=has_legend)
sns.histplot(data=clinvar_zoomed_distances, ax=ax2ins, stat='count', bins=[x/50 for x in range(0,51)], multiple='stack', x='distance', kde=False, hue='label', element='bars', linewidth=0.5, palette=pal, legend=has_legend)
sns.histplot(data=pathogenic_clinvar_zoomed_distances, ax=ax3ins, stat='count', bins=[x/50 for x in range(0,51)], multiple='stack', x='distance', kde=False, hue='label', element='bars', palette=pal, linewidth=0.5, legend=has_legend)
for ax in [ax1ins, ax2ins, ax3ins]:
	ax.spines.top.set_color('grey')
	ax.spines.bottom.set_color('grey')
	ax.spines.left.set_color('grey')
	ax.spines.right.set_color('grey')
	ax.set_xlabel(None)
	ax.set_ylabel(None)
ax1.set_ylabel('GWAS Catalog')
ax2.set_ylabel('Clinvar')
ax3.set_ylabel('Clinvar Pathogenic')
fig.supxlabel('Distance to closest variant (Mbp)')
fig.supylabel('Number of regions')
plt.tight_layout()
fig.savefig('chm13_v1.0_novel_and_nonsyntenic_and_intersection.all.svg', format='svg')