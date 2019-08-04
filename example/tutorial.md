# tutorial

This will walk through the basic analyses required to run phyckle. 

## edge analyses

For edge analyses, we will want to 

- conduct phylogenetic analyses on all the gene trees
- conduct edge calculation analyses
- examine outgroups
- check whether ML trees are consistent with constraint trees
- summarize results

Find the alignment files in the `c1c2_f25_reduced2` folder
`~/phyckle/src/run_raxml.py -r raxml -t 2 -d ../c1c2_f25_reduced2/ -s FNA -p`

Move the ml trees to the `mltrees` folder and the sh-like trees to `shtrees`
`mv RAxML_bestTree.FNA2AA.* mltrees/`
`mv RAxML_fastTreeSH_Support.* shtrees/`
You can remove the rest
`rm RAxML_info* RAxML_log* RAxML_fastTree.* RAxML_parsimonyTree* RAxML_result*`
`rm ../c1c2_f25_reduced2/*.reduced`

Now run the edge investigation :
`~/phyckle/src/edge_investigate_conflicts_given.py -c constraints -o constraints.out -m edge_outs -r raxml -t 2 -d ../c1c2_f25_reduced2/`

This will take a while. 

`~/phyckle/src/postprocess_edge_investigate.py -c constraints.out -d edge_outs`

We can test whether these violate our outgroup requirements. 
`~/phyckle/src/check_outgroup_mono.py -g outgroup.tre -d shtrees -o mono.out -n outgroup`

*todo* Add branch length and taxa filters

We can also test whether the tree that is supported vs our alternatives actually agrees with (or at least doesn't conflict with) the ML tree.
`~/phyckle/src/get_all_ml_constraint_mappings.py -c constraints -d shtrees -s ~/phyckle/src/check_outgroup_mono.py`

Now we can summarize.
`~/phyckle/src/combine_csv.py -f cons_0.csv constraints0.csv constraints1.csv constraints2.csv mono.out -o combined.csv`
`~/phyckle/src/process_combined.py`

The output will look something like this. The output includes the constraint, sum of the DlnL, sum of DlnL for genes with >2lnL differences, number of genes, number of genes with >2lnL differences. The first doesn't include the outgroup and branch length filters. The list of genes after a constraint are the ones that have the highest difference (you can check for errors).

```
constraint sumdiff sum2diff lengenes len2genes
0 265.3902139999973 264.5903229999967 7 6
   238.0962670000008 FNA2AA.5248.removed.realigned.c1c2.f25
   8.36985199999981 FNA2AA.5103.removed.realigned.c1c2.f25
   6.614503999997396 FNA2AA.5165.removed.realigned.c1c2.f25
1 0 0 0 0
2 28.890982000000804 28.890982000000804 2 2
   22.098681999999826 FNA2AA.5233.removed.realigned.c1c2.f25
   6.792300000000978 FNA2AA.5270.removed.realigned.c1c2.f25

outgroup and ml
constraint sumdiff sum2diff lengenes len2genes
0 27.29394699999648 26.494055999995908 6 5
   8.36985199999981 FNA2AA.5103.removed.realigned.c1c2.f25
   6.614503999997396 FNA2AA.5165.removed.realigned.c1c2.f25
   5.1202570000023115 FNA2AA.5161.removed.realigned.c1c2.f25
1 0 0 0 0
2 28.890982000000804 28.890982000000804 2 2
   22.098681999999826 FNA2AA.5233.removed.realigned.c1c2.f25
   6.792300000000978 FNA2AA.5270.removed.realigned.c1c2.f25
```

## combinability tests

For the combinability tests we

- conduct iqtree ML analyses (we use iqtree for several reasons but one is that it already calculates AICc values)
- calculate RF distances between trees (these can be weighted or not)
- conduct the clustering analyses

Let's try the combinability with the same data. First run iqtree like
`~/phyckle/src/run_iqtree.py -i iqtree -t 2 -d ../c1c2_f25_reduced2/`

Create a file with the seq filenames in order
`ls ../c1c2_f25_reduced2/FNA2AA.* > treemap` 

Place them all in a file in the same order
`cat *.treefile > mltrees_comb`

Probably delete extra files
`rm *.bionj *.mldist *.log *.ckp.gz`

Calculate the weighted RF for the trees in the tree file
`bp -rfw -t mltrees_comb > rfw`

Conduct the clustering analysis 
`mv FNA2AA.* ../c1c2_f25_reduced2/`
This will use scaled lengths (-spp) with aic (-a) and stop at weights greater than 4 (-s 4)
`~/phyckle/src/test_clusters.py -d ../c1c2_f25_reduced2/ -m treemap -w rfw -e spp -a -s 4`
Just for testing, we are stopping the test with a graph weight of 4 (none of these support concatenation so there is no need to go longer).

You will get output of (or something similar)
```
../c1c2_f25_reduced2/FNA2AA.5165.removed.realigned.c1c2.f25 mw:0 l:-47260.0841 a:94923.4689 b:96011.8566 f:../c1c2_f25_reduced2/FNA2AA.5165.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5206.removed.realigned.c1c2.f25 mw:0 l:-70087.0566 a:140568.8586 b:141732.4521 f:../c1c2_f25_reduced2/FNA2AA.5206.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5226.removed.realigned.c1c2.f25 mw:0 l:-6541.2451 a:13477.4484 b:14092.2342 f:../c1c2_f25_reduced2/FNA2AA.5226.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5172.removed.realigned.c1c2.f25 mw:0 l:-38584.5889 a:77554.8374 b:78509.3835 f:../c1c2_f25_reduced2/FNA2AA.5172.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5229.removed.realigned.c1c2.f25 mw:0 l:-34109.1368 a:68618.3106 b:69502.3229 f:../c1c2_f25_reduced2/FNA2AA.5229.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5248.removed.realigned.c1c2.f25 mw:0 l:-20275.484 a:40863.332 b:41528.7069 f:../c1c2_f25_reduced2/FNA2AA.5248.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5111.removed.realigned.c1c2.f25 mw:0 l:-28618.025 a:57686.3113 b:58586.927 f:../c1c2_f25_reduced2/FNA2AA.5111.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5233.removed.realigned.c1c2.f25 mw:0 l:-22645.4154 a:45759.5779 b:46654.9887 f:../c1c2_f25_reduced2/FNA2AA.5233.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5103.removed.realigned.c1c2.f25 mw:0 l:-5532.8393 a:11187.9588 b:11430.7508 f:../c1c2_f25_reduced2/FNA2AA.5103.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5270.removed.realigned.c1c2.f25 mw:0 l:-22881.0987 a:46204.3664 b:47000.8514 f:../c1c2_f25_reduced2/FNA2AA.5270.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5276.removed.realigned.c1c2.f25 mw:0 l:-59817.537 a:120065.9584 b:121180.8804 f:../c1c2_f25_reduced2/FNA2AA.5276.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5279.removed.realigned.c1c2.f25 mw:0 l:-19737.5827 a:39884.4351 b:40625.9866 f:../c1c2_f25_reduced2/FNA2AA.5279.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5286.removed.realigned.c1c2.f25 mw:0 l:-30788.1774 a:61943.6592 b:62633.8626 f:../c1c2_f25_reduced2/FNA2AA.5286.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5294.removed.realigned.c1c2.f25 mw:0 l:-16179.7743 a:32763.0896 b:33547.8922 f:../c1c2_f25_reduced2/FNA2AA.5294.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5116.removed.realigned.c1c2.f25 mw:0 l:-12350.8377 a:25274.5624 b:25678.4328 f:../c1c2_f25_reduced2/FNA2AA.5116.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5154.removed.realigned.c1c2.f25 mw:0 l:-22304.8899 a:44916.5826 b:45683.6067 f:../c1c2_f25_reduced2/FNA2AA.5154.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5161.removed.realigned.c1c2.f25 mw:0 l:-46735.1802 a:93913.4703 b:94997.0343 f:../c1c2_f25_reduced2/FNA2AA.5161.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5148.removed.realigned.c1c2.f25 mw:0 l:-19796.6913 a:40025.843 b:40876.9778 f:../c1c2_f25_reduced2/FNA2AA.5148.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5188.removed.realigned.c1c2.f25 mw:0 l:-17230.9897 a:34647.7204 b:35069.9885 f:../c1c2_f25_reduced2/FNA2AA.5188.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5175.removed.realigned.c1c2.f25 mw:0 l:-17081.4354 a:34696.8813 b:35478.8561 f:../c1c2_f25_reduced2/FNA2AA.5175.removed.realigned.c1c2.f25
../c1c2_f25_reduced2/FNA2AA.5170.removed.realigned.c1c2.f25 mw:0 l:-26570.8075 a:53595.5429 b:54567.5488 f:../c1c2_f25_reduced2/FNA2AA.5170.removed.realigned.c1c2.f25
```