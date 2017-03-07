###############################################################################
# Filename: Splice_Site_Plot.py                                               #
# Written by: Lon Blauvelt                                                    #
# Data, stylesheet, and guidance Provided by Prof. Christopher Vollmers       #
#                                                                             #
# Creates a plot of frequencies of nucleotides adjacent to 3'/5' splice sites #
###############################################################################

import matplotlib.pyplot as plt
import sys
import matplotlib.image as mpimg
import numpy as np

# Location:
# /home/lifeisaboutfishtacos/.local/lib/python3.5/site-packages/matplotlib/mpl-data/stylelib
plt.style.use('vollmers')

fig_width      = 4
fig_height     = 1.5
panel_width    = 1.5
panel_height   = 0.75
panel_x_margin = 0.15
panel_y_margin = 0.3
panel_spacing  = 0.025

x_lims = [0,20]
y_lims = [0,2]

plt.rcParams['axes.linewidth'] = 0.4 #set the plot box-width globally

A=mpimg.imread('A.png')
T=mpimg.imread('T.png')
C=mpimg.imread('C.png')
G=mpimg.imread('G.png')

previous_height = 0
base_image = []
x_positions=np.arange(0,20,1)

plt.figure(figsize=(fig_width,fig_height))

###############################################################################
# Read data and prepare it for plotting.                                      #
# Data File: Splice_Sequences.fasta                                           #
###############################################################################

dict5 = {}
dict3 = {}
base_list = ['A', 'T', 'C', 'G']
n_sequences = 0

# Accept filename input from the commandline
# splice_sequence_filename = str(sys.argv[1])
splice_sequence_filename = 'Splice_Sequences.fasta'

for line in open(splice_sequence_filename):
    split_line = line.strip().split('\t')
    splintered_word = split_line[0]
    base_number = 0
    if str(splintered_word[0]) is '>':
        is3or5 = int(splintered_word[1])
    if str(splintered_word[0]) is not '>':
        n_sequences = n_sequences + 1
        if is3or5 is 3:
            for base in splintered_word:
                base_key = str(base_number) + str(base)
                dict3.setdefault(base_key, []).append(1)
                base_number = base_number + 1
        if is3or5 is 5:
            for base in splintered_word:
                base_key = str(base_number) + str(base)
                dict5.setdefault(base_key, []).append(1)
                base_number = base_number + 1

small_sample_correction = (1 / (np.log(2))) * (3 / (2 * n_sequences))

###############################################################################
# Plot #1 (5' Sequences)                                                      #
###############################################################################

# create panel 1
panel1=plt.axes([panel_x_margin, \
                 panel_y_margin, \
                 panel_width/fig_width, \
                 panel_height/fig_height])

# draws a line from (10,0) to (10,2) on panel1
plt.plot([10, 10], [0, 2], 'k-', lw=0.25)

for x in x_positions:
    for base in base_list:
        base_key = str(x) + str(base)
        if base is 'A':
            A_freq = len(dict5[base_key])
        if base is 'T':
            T_freq = len(dict5[base_key])
        if base is 'G':
            G_freq = len(dict5[base_key])
        if base is 'C':
            C_freq = len(dict5[base_key])
    A_relfreq = A_freq / (A_freq + T_freq + G_freq + C_freq)
    T_relfreq = T_freq / (A_freq + T_freq + G_freq + C_freq)
    G_relfreq = G_freq / (A_freq + T_freq + G_freq + C_freq)
    C_relfreq = C_freq / (A_freq + T_freq + G_freq + C_freq)
    H = (A_relfreq * np.log2(A_relfreq)) + \
        (T_relfreq * np.log2(T_relfreq)) + \
        (G_relfreq * np.log2(G_relfreq)) + \
        (C_relfreq * np.log2(C_relfreq))
    information_content = np.log2(4) + H - small_sample_correction
    A_height = A_relfreq * information_content
    T_height = T_relfreq * information_content
    G_height = G_relfreq * information_content
    C_height = C_relfreq * information_content

    heights = sorted([[A_height, A], [T_height, T], [G_height, G], [C_height, C]])
    panel1.imshow(heights[0][1], extent=[x,x+1,0,heights[0][0]],aspect='auto',alpha=1,zorder=1000)
    previous_height+=heights[0][0]
    panel1.imshow(heights[1][1], extent=[x,x+1,previous_height,previous_height+heights[1][0]],aspect='auto',alpha=1,zorder=1000)
    previous_height+=heights[1][0]
    panel1.imshow(heights[2][1], extent=[x,x+1,previous_height,previous_height+heights[2][0]],aspect='auto',alpha=1,zorder=1000)
    previous_height+=heights[2][0]
    panel1.imshow(heights[3][1], extent=[x,x+1,previous_height,previous_height+heights[3][0]],aspect='auto',alpha=1,zorder=1000)
    previous_height=0

# Title of panel1
panel1.text((0.125 + ((panel_width/fig_width))), (0.3 + panel_height), \
             "5'SS", \
             horizontalalignment='center', \
             verticalalignment='center', \
             transform = panel1.transAxes)

# Set Axes Limits
panel1.set_xlim(x_lims)
panel1.set_ylim(y_lims)

# Set Axes Ticks
panel1.set_xticks([0,5,10,15,20])
panel1.set_yticks([0.0,0.5,1.0,1.5,2.0])

# Set Axes Tick Labels
panel1.set_xticklabels([-10,-5,0,5,10])

# Set Axes Tick Thickness
panel1.xaxis.set_tick_params(width=0.4)
panel1.yaxis.set_tick_params(width=0.4)

# Set Axes Labels
panel1.set_ylabel("Bits")
panel1.set_xlabel("Distance to\nSplice Site")

# Turn On or Off Plot Axes
panel1.tick_params(axis='both', which='both', \
                   bottom='on', labelbottom='on', \
                   left='on', labelleft='on', \
                   right='off', labelright='off', \
                   top='off', labeltop='off', \
                   length=1)

###############################################################################
# Plot #2 (3' Sequences)                                                      #
###############################################################################

panel2=plt.axes([panel_x_margin + panel_width/fig_width + panel_spacing, \
                 panel_y_margin, \
                 panel_width/fig_width, \
                 panel_height/fig_height])

# draws a line from (10,0) to (10,2) on panel2
plt.plot([10, 10], [0, 2], 'k-', lw=0.25)

for x in x_positions:
    for base in base_list:
        base_key = str(x) + str(base)
        if base is 'A':
            A_freq = len(dict3[base_key])
        if base is 'T':
            T_freq = len(dict3[base_key])
        if base is 'G':
            G_freq = len(dict3[base_key])
        if base is 'C':
            C_freq = len(dict3[base_key])

    A_relfreq = A_freq / (A_freq + T_freq + G_freq + C_freq)
    T_relfreq = T_freq / (A_freq + T_freq + G_freq + C_freq)
    G_relfreq = G_freq / (A_freq + T_freq + G_freq + C_freq)
    C_relfreq = C_freq / (A_freq + T_freq + G_freq + C_freq)
    H = (A_relfreq * np.log2(A_relfreq)) + \
        (T_relfreq * np.log2(T_relfreq)) + \
        (G_relfreq * np.log2(G_relfreq)) + \
        (C_relfreq * np.log2(C_relfreq))
    information_content = np.log2(4) + H - small_sample_correction
    A_height = A_relfreq * information_content
    T_height = T_relfreq * information_content
    G_height = G_relfreq * information_content
    C_height = C_relfreq * information_content

    heights = sorted([[A_height, A], [T_height, T], [G_height, G], [C_height, C]])
    panel2.imshow(heights[0][1], extent=[x,x+1,0,heights[0][0]],aspect='auto',alpha=1,zorder=1000)
    previous_height+=heights[0][0]
    panel2.imshow(heights[1][1], extent=[x,x+1,previous_height,previous_height+heights[1][0]],aspect='auto',alpha=1,zorder=1000)
    previous_height+=heights[1][0]
    panel2.imshow(heights[2][1], extent=[x,x+1,previous_height,previous_height+heights[2][0]],aspect='auto',alpha=1,zorder=1000)
    previous_height+=heights[2][0]
    panel2.imshow(heights[3][1], extent=[x,x+1,previous_height,previous_height+heights[3][0]],aspect='auto',alpha=1,zorder=1000)
    previous_height=0

# Title of panel2
panel2.text((0.125 + ((panel_width/fig_width))), (0.3 + panel_height), \
             "3'SS", \
             horizontalalignment='center', \
             verticalalignment='center', \
             transform = panel2.transAxes)

# Set Axes Limits
panel2.set_xlim(x_lims)
panel2.set_ylim(y_lims)

# Set Axes Ticks
panel2.set_xticks([0,5,10,15,20])
panel2.set_yticks([0.0,0.5,1.0,1.5,2.0])

# Set Axes Tick Labels
panel2.set_xticklabels([-10,-5,0,5,10])

# Set Axes Tick Thickness
panel2.xaxis.set_tick_params(width=0.4)
panel2.yaxis.set_tick_params(width=0.4)

# Set Axes Label
panel2.set_xlabel("Distance to\nSplice Site")

# Turn On or Off Plot Axes
panel2.tick_params(axis='both', which='both', \
                   bottom='on', labelbottom='on', \
                   left='off', labelleft='off', \
                   right='off', labelright='off', \
                   top='off', labeltop='off', \
                   length=1)

plt.savefig('Plot_Splice_Frequency.pdf')
