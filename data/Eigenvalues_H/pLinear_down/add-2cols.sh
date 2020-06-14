#!/usr/bin/bash

for f in *.dat
do
  # backup...
  cp ${f} ${f/.dat/.dtt}
  # remove DOS...
  sed -i -e 's@@@' ${f}
  head -n 2 ${f} > ${f/.dat/.ttt}
  grep -v '^#' ${f} | awk 'BEGIN { OFMT="%.6f" } { print $0, ($5/$2), ($6/$2) ; }' | sed -e 's@ @	@g' >> ${f/.dat/.ttt}
  mv ${f/.dat/.ttt} ${f}
done

# vim: noet list

