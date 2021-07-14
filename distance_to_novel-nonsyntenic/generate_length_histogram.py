import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import numpy as np


plt.rcParams['svg.fonttype'] = 'none'

#=====================#
# Define region files #
#=====================#

# These are the raw novel and non-syntenic region files
# Regions may overlap between files

novel_regions_file = 'chm13_v1.0_uncoveredByGRCh38WinnowmapAlignments.bed'
nonsyntenic_regions_file = 'chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp.bed'


#======================#
# Load in region sizes #
#======================#

region_sizes_by_type = {}
zoomed_sizes_by_type = {}

for line in open(novel_regions_file):
	if line.startswith('#'):
		continue
	fields = line.rstrip('\n').split('\t')
	size = int(fields[2]) - int(fields[1])
	region_sizes_by_type.setdefault('NOVEL', [])
	region_sizes_by_type['NOVEL'].append(size)


for line in open(nonsyntenic_regions_file):
	if line.startswith('#'):
		continue
	fields = line.rstrip('\n').split('\t')
	size = int(fields[2]) - int(fields[1])
	region_sizes_by_type.setdefault('NONSYNTENIC', [])
	region_sizes_by_type['NONSYNTENIC'].append(size)


# This converts regions lengths to a log base 10 scale
for region in region_sizes_by_type:
	region_sizes_by_type[region] = np.log10(region_sizes_by_type[region])


#==================#
# Generate Figures #
#==================#

sns.set_theme()

bins = [x/10 for x in range(0,81,2)]

fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(10,5))
sns.histplot(region_sizes_by_type['NOVEL'], ax=ax1, stat='count', kde=False, bins=bins, linewidth=0.5, element='bars')
sns.histplot(region_sizes_by_type['NONSYNTENIC'], ax=ax2, stat='count', kde=False, bins=bins, linewidth=0.5, element='bars')
ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
fig.supxlabel('Size of Region (log10[bp])')
ax1.set_ylabel('Novel')
ax2.set_ylabel('Non-syntenic')
fig.supylabel('Number of Regions')
plt.tight_layout()
fig.savefig('region_sizes_hist.svg', dpi=400)
