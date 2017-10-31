#!/usr/bin/env python3
'''
svgutils_BESTD.py
=================

:Author: |author_name|
:Release: |version|
:Date: |today|


Purpose
=======

|description|

Make a figure layout using the python package svgutils.
These scripts are hard coded for S1, S2 and S7 plots of the BESTD paper.


See examples and information from svgutils:

    - https://neuroscience.telenczuk.pl/?p=331
    - http://svgutils.readthedocs.io/en/latest/tutorials/composing_multipanel_figures.html
    - https://svgutils.readthedocs.io/en/latest/compose.html
    - http://cairosvg.org/
    - https://inkscape.org/en/


Input:

    Requires svg files as input.


Output:

    Outputs svg and pdf formats of a multi-panel plot.

Requires:

    Python packages and Inkscape


Documentation
=============

    For more information see:

        |url|

'''
##############
# Get all the modules needed
# System:
import os
import sys
import glob

# Options and help:
import docopt

# Get svgutils and PDF converter if not using Inkscape:
from svgutils.compose import *
import cairosvg

# Get additional packages needed:
import string
##############

##############
plotA = 'S2A_plot_PCA_raw.svg'
plotB = 'S2B_density_plot_normalised_VSN.svg'
plotC = 'S2C_SD_vs_means_normalised_GEx.svg'
figure_name = 'S2_GEx_QC'

pos1 = 5
pos2 = 5
moves = [(20, 20), (380, 20), (20, 380), (380, 380)]
scales = [0.65, 0.75, 0.8, 0.7]
size = 11
weight = 'bold'

file_format_in = 'svg'
file_format_out = 'pdf'
layout_name_1 = str('{}.{}'.format(figure_name,
                                   file_format_in))
layout_name_2 = str('{}.{}'.format(figure_name,
                                   file_format_out))
##############


##############
fig_layout =  Figure("21cm", "19.7cm",
                     # Plot A, then B, etc. use scale and move for each panel
                     # individually:
                     Panel(SVG(plotA).scale(scales[0]),#.move(moves[0][0], moves[0][1]),
                           # scale only the plot, not the text
                           Text("A", pos1, pos2, size = size, weight = weight),
                           # Place Text() after SVG(), otherwise it doesn't plot
                           ).move(moves[0][0], moves[0][1]),
                           # move the whole panel
                     Panel(SVG(plotB).scale(scales[1]),#.move(moves[1][0], moves[1][1]),
                           Text("B", pos1, pos2, size = size, weight = weight),
                           ).move(moves[1][0], moves[1][1]),
                     Panel(SVG(plotC).scale(scales[2]).move(10, 10),
                           Text("C", pos1, pos2, size = size, weight = weight),
                           ).move(moves[2][0], moves[2][1]),
                     #Panel(SVG(plotD).scale(scales[3])..move(moves[3][0], moves[3][1]),
                     #      Text("D", pos1, pos2, size = size, weight = weight),
                     #      ).move(moves[3][0], moves[3][1])
                     )

# Save the Figure:
fig_layout.save(layout_name_1)


os.system('''inkscape --without-gui \
                              --export-area-drawing \
                              --export-margin=2 \
                              --file={} \
                              --export-pdf={}'''.format(layout_name_1,
                                                        layout_name_2)
         )
##############


##############
# Finish and exit with docopt arguments:
#if __name__ == '__main__':
#    sys.exit(main())
############
