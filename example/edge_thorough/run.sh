tar -xzf seqs.tar.gz
SEQD=seqs
PHYC=~/phyckle/src
PHYC=~/Dropbox/programming/phyckle/src
$PHYC/run_raxml.py -r raxml -t 2 -d $SEQD/ -s fa -p
mkdir mltrees shtrees
mv RAxML_bestTree.* mltrees/
mv RAxML_fastTreeSH_Support.* shtrees/
rm RAxML_info* RAxML_log* RAxML_fastTree.* RAxML_parsimonyTree* RAxML_result*
rm $SEQD/*.reduced

#first one
mkdir con_first
cd con_first
$PHYC/edge_investigate_conflicts_given.py -c ../constraints_first -o constraints.out -m edge_outs -r raxml -t 2 -d ../$SEQD/
rm RAxML_log.g* RAxML_result.g*
$PHYC/postprocess_edge_investigate.py -c constraints.out -d edge_outs
$PHYC/get_all_ml_constraint_mappings.py -c ../constraints_first -d ../shtrees -s $PHYC/check_outgroup_mono.py
$PHYC/combine_csv.py -f cons_0.csv constraints0.csv constraints1.csv constraints2.csv -o combined.csv
$PHYC/process_combined.py -f combined.csv
cd ../

#second one
mkdir con_second
cd con_second
$PHYC/edge_investigate_conflicts_given.py -c ../constraints_second -o constraints.out -m edge_outs -r raxml -t 2 -d ../$SEQD/
rm RAxML_log.g* RAxML_result.g*
$PHYC/postprocess_edge_investigate.py -c constraints.out -d edge_outs
$PHYC/get_all_ml_constraint_mappings.py -c ../constraints_second -d ../shtrees -s $PHYC/check_outgroup_mono.py
$PHYC/combine_csv.py -f cons_0.csv constraints0.csv constraints1.csv constraints2.csv -o combined.csv
$PHYC/process_combined.py -f combined.csv
cd ../
