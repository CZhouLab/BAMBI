#! /bin/bash

head -1 clusters.w_lnc_biomarkers   > cluster1.row
row2col_space.pl cluster1.row > cluster1.col

for i in {2..2} 
do
	head -$i clusters.w_lnc_biomarkers | tail -1 > cluster$i.row
	row2col_space.pl cluster$i.row > cluster$i.col
done
