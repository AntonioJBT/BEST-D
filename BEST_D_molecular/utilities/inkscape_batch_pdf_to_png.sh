#!/usr/bin/env bash

# run as e.g. bash inkscape_batch_pdf_to_png.sh PCA_SNVs_all_samples.pdf
# For a loop do e.g.:
# for i in *pdf; do echo $i; done
# for i in *pdf; do ./inkscape_batch_pdf_to_png.sh $i; done

# Set bash script options:
# https://kvz.io/blog/2013/11/21/bash-best-practices/
set -o errexit
set -o pipefail
set -o nounset


INKSCAPE_PATH=inkscape
SOURCE_FILE=$1
DEST_PNG=${SOURCE_FILE}
DPI=600

${INKSCAPE_PATH} --without-gui -f ${SOURCE_FILE} --export-png ${DEST_PNG} --export-area-drawing --export-dpi=${DPI}.png
