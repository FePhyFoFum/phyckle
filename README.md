# phyckle

In the `src` folder are several scripts to help conduct phylogenomic analyses.

### conducting combinability analyses

At some point for this, you will need to use `iqtree` and `bp` (from [here](https://github.com/FePhyFoFum/gophy)).

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

At some point these analyses will require `RAxML`. 

There are several ways that these analyses can be conducted. Once means to compare alternatives is to do the following:

- Create a file that has the alternative bipartitions listed. There is an example in `example/test.bp` with the focal bipartition on a line and then alternatives listed below preceded by a `-`.
- Run `edge_investigate_conflicts_given.py` like `python phyckle/src/edge_investigate_conflicts_given.py test.bp test.out seqs/ temp/` where `temp` is an outdir, `seqs` only has the sequence files (nothing else), `test.bp` is the set of bipartitions and `test.out` is an outfile. 
- Then run `postprocess_edge_investigate.py` link `python phyckle/src/postprocess_edge_investigate.py test.out temp/` with `test.out` from above and `temp` from above. 
- This will result in a file called `cons_0.csv` that has each line with `gene,cons_0,cons_0_conf_0,cons_0_conf_1,bestone,best,secondbest,diffbestsecondbest` listing the `gene`, likelihoods for each edge (`cons_0` is the main and `cons_0_conf_0,cons_0_conf_1` are the two alternatives from an example). `best` is the highest likelihood with the differences listed.

There are several other analyses that can be conducted but this is the basic set of analyses. 

### some comments

The runs above switch between `iqtree` and `RAxML`. This is primarily because `iqtree` has several branch length options when concatenating and `RAxML` is more efficient when conducting constrained analyses. 