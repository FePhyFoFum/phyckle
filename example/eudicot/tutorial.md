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