#!/bin/bash

# run dot on all dot files
# The DOT language is a plain text graph language:
# http://www.graphviz.org/content/dot-language
# If you have doxygen installed you probably have 
# already installed dot.
# see forward_operator.gv for a simple example.

for f in *.gv

do

  dot $f -Tsvg -o $f'.svg'
  #dot $f -Tpdf -o $f'.pdf'

done


