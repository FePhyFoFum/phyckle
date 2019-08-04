mkdir seqs
for i in {1..10}
do
    pxtscale -t t1 -r 0.75 | pxseqgen > g$i.fa
done
for i in {10..15}
do
    pxtscale -t t2 -r 0.75 | pxseqgen > g$i.fa
done
for i in {15..20}
do
    pxtscale -t t3 -r 0.75 | pxseqgen > g$i.fa
done
mv g[0-9]*.fa seqs/