t_coffee -reg -reg_method famsa_msa \
     -seq ${seqs} \
     -reg_tree ${guide_tree} \
     -child_tree ${slave_tree} \
     -reg_nseq ${bucket_size} \
     -reg_homoplasy \
     -outfile ${id}.slave_align.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.aln