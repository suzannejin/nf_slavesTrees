t_coffee -other_pg seq_reformat -in $seqs -in2 $guide_tree -action +regtrim $nseqid > ${id}.informative${nseqid}.with.${tree_method}.tree.fa

t_coffee -other_pg seq_reformat -in ${id}.informative${nseqid}.with.${tree_method}.tree.fa -output sim_idscore > ${id}.informative${nseqid}.with.${tree_method}.tree.sim_idscore
awk '\$1=="AVG"&&\$2=="AVG"{print \$4}' ${id}.informative${nseqid}.with.${tree_method}.tree.sim_idscore > ${id}.informative${nseqid}.with.${tree_method}.tree.seqidscore
