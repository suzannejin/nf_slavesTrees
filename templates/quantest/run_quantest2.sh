python3 $baseDir/bin/get_refANDinformative_seqs.py \
    --msa ${msa} \
    --seq ${seqs} \
    --names ${refnames} \
    --tree ${guide_tree} \
    --n ${nquan} \
    >> ${id}.slave_align.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.informative${nquan}.with.${tree_method}.tree.aln

quantest2 ${id}.slave_align.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.informative${nquan}.with.${tree_method}.tree.aln ${ss} 

awk 'NR==4{print}' ${id}.slave_align.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.informative${nquan}.with.${tree_method}.tree.quantest2_log \
    | cut -f2 \
    > ${id}.slave_align.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.informative${nquan}.with.${tree_method}.tree.quantest2

