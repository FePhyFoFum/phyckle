# phyckle

In the `src` folder are several scripts to help conduct phylogenomic analyses.

### conducting combinability analyses

There are two scripts for conducting these analyses currently: `test_clusters.py` and `test_clusters_multi_measure.py`. The `test_clusters.py` will conduct analyses by using the RFW, constructing a graph, and proceeding through the graph. `test_clusters_multi_measure.py` has additional considerations like penalizing based on poor overlap of taxa (using `rfp` from `bp`). 

- Conduct iqtree analyses on each gene. You can do this with `run_iqtree.py`.
- Create a file with the seq filenames in order (`ls *.fa > treemap`) and place them all in a file in the same order (`cat *.treefile > mltrees`)
- Probably delete extra files (`rm *.bionj *.mldist *.log *.ckp.gz`) and move iqtree files and treefiles to the gene directory.
- Calculate the weighted RF for the trees in the tree file (`bp -rfw -t mltres > rfw`).
- Conduct the clustering analysis `python phyckle/src/test_clusters.py -d seqs/ -m treemap -w rfw -e spp -a`

_Example with pxbdsim and pxseqgen_. You can use phyx to generate data that you can use to conduct an example analysis. Here are the commands:
```
pxbdsim -e 10 | pxtscale -r 1 > test.tre
pxseqgen -t test.tre -l 500 -o test1.fa
pxseqgen -t test.tre -l 500 -o test2.fa
pxseqgen -t test.tre -l 500 -o test3.fa
pxseqgen -t test.tre -l 500 -o test4.fa
pxseqgen -t test.tre -l 500 -o test5.fa
mkdir seqs/
mv *.fa seqs/
python phyckle/src/run_iqtree.py seqs/ test
rm *.bionj *.mldist *.log *.ckp.gz
mv *.treefile *.iqtree seqs/
cd seqs/
ls test*.fa > ../treemap
cat test*.treefile > ../mltrees
cd ../
bp -rfw -t mltrees > rfw
python phyckle/src/test_clusters.py -d seqs/ -m treemap -w rfw -e spp -a
```

### conducting alternative relationship analyses