import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import numpy as np

plt.rcParams['svg.fonttype'] = 'none'

#=====================#
# Define region files #
#=====================#

# This is a bed file containing the non-overlapping, annotated novel and non-syntenic regions

novel_regions_file = 'chm13_v1.0_novel_and_nonsyntenic_and_intersection.all.sorted.bed'


#========================#
# Load in region lengths #
#========================#

region_sizes_by_type = {}
zoomed_sizes_by_type = {}

for line in open(novel_regions_file):
	if line.startswith('#'):
		continue
	fields = line.rstrip('\n').split('\t')
	size = int(fields[2]) - int(fields[1])
	region_type = fields[3]
	region_sizes_by_type.setdefault(region_type, [])
	region_sizes_by_type[region_type].append(size)

# This converts regions lengths to a log base 10 scale
for region_type in region_sizes_by_type:
	region_sizes_by_type[region_type] = np.log10(region_sizes_by_type[region_type])


#==================#
# Generate Figures #
#==================#

sns.set_theme()


bins = [x/10 for x in range(0,81,2)]

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True, figsize=(10,7))
sns.histplot(region_sizes_by_type['NOVEL'], ax=ax1, stat='count', kde=False, bins=bins, linewidth=0.5, element='bars')
sns.histplot(region_sizes_by_type['NONSYNTENIC'], ax=ax2, stat='count', kde=False, bins=bins, linewidth=0.5, element='bars')
sns.histplot(region_sizes_by_type['BOTH'], ax=ax3, stat='count', kde=False, bins=bins, linewidth=0.5, element='bars')
fig.supxlabel('Size of Region (log10[bp])')
ax1.set_ylabel('Novel Only')
ax2.set_ylabel('Non-syntenic Only')
ax3.set_ylabel('Both')
fig.supylabel('Number of Regions')
plt.tight_layout()
fig.savefig('no_overlap_region_sizes_hist.svg', dpi=400)