#!/usr/bin/env python

#Program to plot OTU abundances for each sample in an otu table, and to plot
#OTU richness curves.
#Author: Dr. Benjamin C. Smith, 2011

#Usage: sample_otu_plots_v0.1.py pooled_table.txt out_dir/


import re, sys, os
from numpy import array, loadtxt, vstack, hstack, arange
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages


def sample_otus(table, column=1, threshold=0):
    """Extracts list of OTUs in sample number given by col, that have more
    than the number of reads set by threshold.
    """
    
    otus = []
    for i, otu_ct in enumerate(table[:, column]):
        if float(otu_ct) >= threshold:
            label = table[i, -1]#.split("_")[0]
            otus.append([label, int(otu_ct)])
    return len(otus), otus


file = str(sys.argv[1])
print file
otu_table = loadtxt(file, dtype = 'object')

data = open(file, 'r')
for n, line in enumerate(data):
    if n == 1:
        headers = line.split("\t")[1:-1]

step = 10.
threshold_proportion = 0.005
x = arange(10., 1001., step)
fig_nums = []
pie_nums = []
for i, name in enumerate(headers[:]):
    otu_counts = []
    #get number of counts at a range of count thresholds
    for j, lim in enumerate(x):
        sample_otus_out = sample_otus(otu_table, column=i+1, threshold=lim)
        counts, otus = sample_otus_out
        otu_counts.append(counts)
        #pick a count threshold at which to generate pie charts
        if len(otus) > 0:
            all_otus = vstack( array( otus, dtype=object) )
            max_otu_ct = all_otus[:,1].max()
        #plot_threshold = (step/2) + (max_otu_ct * threshold_proportion)
        plot_threshold = 50
        if plot_threshold > x[j] and plot_threshold <= x[j+1] :
            otus_4_pie = otus
    #plot otu richness figures with "ncurves" curves to a figure
    ncurves = 5
    fnum = int(floor(i / ncurves) + 1)
    fig_nums.append(fnum)
    if array( otu_counts ).max() != 0:
        figure(fnum)
        subplot(111)
        semilogy(x, otu_counts, label=name)
        legend(bbox_to_anchor=(0.50, .98), loc=2, borderaxespad=0.)
        title("OTU Richness Plot")
        ylabel("Number of OTUs")
        xlabel("Acceptance threshold (counts)")
    #plot pie charts
    pnum = int((len(headers)/ncurves)+i+10)
    pie_nums.append(pnum)
    figure(pnum, figsize=(8,8))
    ax = axes([0.1, 0.1, 0.8, 0.8])
    subplot(2,2,1)
    otus_4_pie.sort(lambda x, y: cmp(x[1],y[1]), reverse=True)
    max_otus = len(otus_4_pie)
    if len(otus_4_pie) > max_otus:
        pie_otus = otus_4_pie[0:max_otus]
        others = 0
        for num in otus_4_pie[max_otus:]:
            others += num[1]
        pie_otus.append(["Others", others])
    else:
        pie_otus = otus_4_pie
    fracs = []
    for line in pie_otus:
        fracs.append(line[1])
    pie(fracs, autopct='%1.1f%%', pctdistance=1.25)
    title(headers[i])
    pie_labels = []
    for line in pie_otus:
        pie_labels.append(line[0] + " (cts=" + str(line[1]) + ")")
    legend(pie_labels, loc=3, bbox_to_anchor=(0., -.9))



#save figures and pie charts
for i in fig_nums:
    figure(i)
    savename = sys.argv[2] + "otu_richness_" + str(i) + ".pdf"
    savefig(savename, format='pdf')

for i, num in enumerate(pie_nums):
    figure(num)
    savename = sys.argv[2] + headers[i] + " Pie.pdf"
    savefig(savename, format='pdf')



# create figure
# figwidth = 8.5    # inches
# figheight = 6.0   # inches
# figure(1, figsize=(figwidth, figheight))
# rcParams['font.size'] = 12.0
# rcParams['axes.titlesize'] = 16.0
# rcParams['xtick.labelsize'] = 12.0
# rcParams['legend.fontsize'] = 10.0
# 
# figure.plot()
# 
# 
# show()



    
# for i in range(len(d)): 
#     m = re.findall('Qv[0-9]+rpt_[0-9]+|Qv[0-9]+_[0-9]+', d[i])
#     c.append(len(m))
# 
# 
# otu_id = np.arange(0, len(d))
# pos = otu_id + .5
# val_ticks = np.arange(0, len(d), 5)
# pos_ticks = val_ticks + .5
# 
# plt.barh(pos, c, align='center')
# plt.xlabel('Number of reads')
# plt.yticks(pos_ticks, val_ticks)
# plt.ylabel('OTU')
# plt.axis([0, 12500, 0, len(d)])
# plt.title('Number of reads assigned to each OTU')
# plt.show()
