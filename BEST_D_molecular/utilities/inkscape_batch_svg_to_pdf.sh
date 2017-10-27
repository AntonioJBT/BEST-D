#!/usr/bin/env bash
  
# run as e.g. bash inkscape_batch_svg_to_pdf.sh PCA_SNVs_all_samples.svg
# For a loop do e.g.:
# for i in *svg; do echo $i; done
# for i in *svg; do ./inkscape_batch_svg_to_pdf.sh $i; done

# Set bash script options:
# https://kvz.io/blog/2013/11/21/bash-best-practices/
set -o errexit
set -o pipefail
set -o nounset


INKSCAPE_PATH=inkscape
SOURCE_FILE=$1
DEST=${SOURCE_FILE}

${INKSCAPE_PATH} --without-gui --file=${SOURCE_FILE} --export-pdf=${DEST}.pdf --export-area-drawing --export-margin=2

rename -s .svg.pdf .pdf *.svg.pdf
