#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'XXXXXX'.
 *
 *   XXXXXX is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   XXXXXX is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with XXXXXX.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main XXX pipeline script
 *
 * @authors
 * Edgar Garriga

 */

/*
 * defaults parameter definitions
 */

// input sequences to align in fasta format
//params.seqs ="${baseDir}/test/*.fa"
params.seqs ="/users/cn/egarriga/datasets/homfam/combinedSeqs/*.fa"

// input reference sequences aligned in
//params.refs ="${baseDir}/test/*.ref"
params.refs ="/users/cn/egarriga/datasets/homfam/refs/*.ref"

// input guide trees in Newick format. Or `false` to generate trees
//params.trees ="/home/edgar/CBCRG/nf_homoplasty/results_trees"

// which tree methods to run if `trees` == `false`
params.tree_method = "CLUSTALO,MAFFT_PARTTREE"
//"famsaUpgma,famsaSL,famsaParttreeSL,famsaParttreeUpgma"     //FAMSA,CLUSTALO,MAFFT_PARTTREE,dpparttreednd0
//codnd,dpparttreednd0,dpparttreednd1,dpparttreednd2,dpparttreednd2size,fastaparttreednd,fftns1dnd,fftns1dndmem,fftns2dnd,fftns2dndmem,mafftdnd,parttreednd0,parttreednd1,parttreednd2,parttreednd2size

// which alignment methods to run
params.align_method = "CLUSTALO"      //"CLUSTALO,MAFFT-FFTNS1,MAFFT-SPARSECORE,UPP,MAFFT-GINSI"

// bucket sizes for regressive algorithm
params.buckets= '1000'

//run reg with slave trees
params.slave_align = true
params.slave_tree_method = "mbed,parttree,famsadnd"

// evaluate alignments ?
params.evaluate = true
params.homoplasy = false
params.metrics = false
params.esl = false

// output directory
params.outdir = "$baseDir/results"

log.info """\
         F  A  M  S  A    P  i  p  e  l  i  n  e  ~  version 0.1"
         ======================================="
         Input sequences (FASTA)                        : ${params.seqs}
         Input references (Aligned FASTA)               : ${params.refs}
         Input trees (NEWICK)                           : ${params.trees}
         Alignment methods                              : ${params.align_method}
         Tree methods                                   : ${params.tree_method}
         Generate Slave tree alignments                 : ${params.slave_align}
         Slave Tree methods                             : ${params.slave_tree_method}
         Bucket size                                    : ${params.buckets}
         Perform evaluation? Requires reference         : ${params.evaluate}
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()


// Channels containing sequences
if ( params.seqs ) {
  Channel
  .fromPath(params.seqs)
  .map { item -> [ item.baseName, item] }
  .into { seqsCh; seqs2 }
}

if ( params.refs ) {
  Channel
  .fromPath(params.refs)
  .map { item -> [ item.baseName, item] }
  .set { refs }
}

// Channels for user provided trees or empty channel if trees are to be generated [OPTIONAL]
if ( params.trees ) {
  Channel
    .fromPath(params.trees)
    .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item] }
    .set { trees }
}
else { 
  Channel
    .empty()
    .set { trees }
}

align_methods = params.align_method
slave_trees_methods = params.slave_tree_method
tree_methods = params.tree_method

/*
 * GENERATE GUIDE TREES USING MEHTODS DEFINED WITH "--tree_method"
 *
 * NOTE: THIS IS ONLY IF GUIDE TREES ARE NOT PROVIDED BY THE USER
 * BY USING THE `--trees` PARAMETER
 */

process generate_trees {
    tag "${id}.${tree_method}"
    publishDir "${params.outdir}/guide_trees", mode: 'copy', overwrite: true
   
    input:
    set val(id), \
         file(seqs) \
         from seqsCh
    each tree_method from tree_methods.tokenize(',')

   output:
     set val(id), \
       val(tree_method), \
       file("${id}.${tree_method}.dnd") \
       into treesGenerated

   when:
     !params.trees

   script:
    template "tree/generate_tree_${tree_method}.sh"

}

treesGenerated
  .mix ( trees )
  .combine ( seqs2, by:0 )
  .into {seqsAndTreesForRegressiveAlignment; seqsAndTreesForSlaveAlignment }


process slave_alignment {
    tag "${id}"
    publishDir "${params.outdir}/alignments", pattern: '*.aln', mode: 'copy', overwrite: true

    input:
      set val(id), \
        val(tree_method), \
        file(guide_tree), \
        file(seqs) \
        from seqsAndTreesForSlaveAlignment

      each slave_tree from slave_trees_methods.tokenize(',')   

      each align_method from align_methods.tokenize(',')   

      each bucket_size from params.buckets.tokenize(',')

    when:
      params.slave_align

    output:
      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val("slave_align"), \
        val(bucket_size), \
        val("${slave_tree}"), \
        file("*.aln") \
        into slaveOut
      
      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val(bucket_size), \
        val("${slave_tree}"), \
        file("${id}.homoplasy") \
        into homoSlave

      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val("slave_align"), \
        val(bucket_size), \
        val("${slave_tree}"), \
        file(".command.trace") \
        into metricsSlave

   script:
    template "slave_reg/slave_reg_${align_method}.sh"
}

process metrics{
    tag "${id}"
    publishDir "${params.outdir}/metrics", mode: 'copy', overwrite: true

    input:
    set val(id), \
      val(align_method), \
      val(tree_method), \
      val(mode), \
      val(bucket_size), \
      val(slave_tree), \
      val(metricsFile) \
      from metricsSlave

    when:
      params.metrics

    output:
    set file("${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.metrics"), \
      file("*.realtime"), \
      file("*.rss"), \
      file("*.peakRss"), \
      file("*.vmem"), \
      file("*.peakVmem"), \
      file("*.metrics") \
        into metricsOut

    script:
    """    
    ## realtime > Task execution time i.e. delta between completion and start timestamp i.e. compute wall-time
    awk -F = '{ if (\$1=="realtime") print \$2}' ${metricsFile} > ${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.realtime

    ## rss > Real memory (resident set) size of the process
    awk -F = '{ if (\$1=="rss") print \$2}' ${metricsFile}> ${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.rss

    ## peakRss > Peak of real memory
    awk -F = '{ if (\$1=="peak_rss") print \$2}' ${metricsFile} > ${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.peakRss

    ## vmem > Virtual memory size of the process
    awk -F = '{ if (\$1=="vmem") print \$2}' ${metricsFile} > ${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.vmem

    ## peakVmem > Peak of virtual memory
    awk -F = '{ if (\$1=="peak_vmem") print \$2}' ${metricsFile} > ${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.peakVmem

    mv ${metricsFile} ${id}.${mode}.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.metrics
    """
}

process homoplasy{
    tag "${id}"
    publishDir "${params.outdir}/homoplasy", mode: 'copy', overwrite: true

    input:
    set val(id), \
      val(align_method), \
      val(tree_method), \
      val(bucket_size), \
      val(slave_tree), \
      file(homoplasy) \
      from homoSlave

    when:
      params.homoplasy

    output:
    set file("*.homo"), \
        file("*.w_homo"), \
        file("*.w_homo2"), \
        file("*.len"), \
        file("*.ngap"), \
        file("*.ngap2") \
        into homoplasyOut

    script:
    """    
    ## homo
    awk -F : '{ if (\$1=="HOMOPLASY") print \$2}' ${homoplasy} > ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.homo
    ## w_homo
    awk -F : '{ if (\$1=="WEIGHTED_HOMOPLASY") print \$2}' ${id}.homoplasy > ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.w_homo
    ## w_homo2
    awk -F : '{ if (\$1=="WEIGHTED_HOMOPLASY2") print \$2}' ${id}.homoplasy > ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.w_homo2
    ## len
    awk -F : '{ if (\$1=="LEN") print \$2}' ${id}.homoplasy > ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.len
    ## ngap
    awk -F : '{ if (\$1=="NGAP") print \$2}' ${id}.homoplasy > ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.ngap
    ## ngap2
    awk -F : '{ if (\$1=="NGAP2") print \$2}' ${id}.homoplasy > ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.ngap2 
    """
}

refs
  .cross (slaveOut )
  .map { it -> [it[0][0], it[1][1], it[1][2], it[1][3], it[1][4], it[1][5], it[1][6],it[0][1]] }
  .into { toEsl; toEvaluate}

process esl{
    tag "${id}.${align_method}.${tree_method}.${align_type}.${bucket_size}"
    publishDir "${params.outdir}/esl", mode: 'copy', overwrite: true
    label 'process_low'
    container 'edgano/hmmer'

    input:
      set val(id), \
          val(align_method), \
          val(tree_method), \
          val(align_type), \
          val(bucket_size), \
          val(slave_tree), \
          file(test_alignment), \
          file(ref_alignment) \
          from toEsl
    when:
      params.esl
    output:
      set file("*.easel_INFO"),file("*.avgLen"),file("*.avgId") into eslOut
      
     shell:
     '''
     esl-alistat !{test_alignment} > !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.slave.!{slave_tree}.easel_INFO
     awk -F : '{ if (\$1=="Average length") print \$2}' !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.slave.!{slave_tree}.easel_INFO | sed 's/ //g' > !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.slave.!{slave_tree}.avgLen 
     awk -F : '{ if (\$1=="Average identity") print substr(\$2, 1, length(\$2)-1)}' !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.slave.!{slave_tree}.easel_INFO | sed 's/ //g' > !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.slave.!{slave_tree}.avgId 

     ## awk 'NR > 8 && $1 !~/\\// { sum+= $3 } END {print "SUM: "sum"\\nAVG: "sum/(NR-9)}' !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.slave.!{slave_tree}.easel_INFO > !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.slave.!{slave_tree}.easel_AVG
     ## the first && is to skip first lines and the last one. The AVG is done -8 all the time execpt for the END print to "erase" the last "//" too.
     '''
}

process evaluation {
    tag "${id}.${align_method}.${tree_method}.${align_type}.${bucket_size}"
    publishDir "${params.outdir}/individual_scores", mode: 'copy', overwrite: true
    //container 'edgano/homoplasy:latest'
    label 'process_low'

    input:
      set val(id), \
          val(align_method), \
          val(tree_method), \
          val(align_type), \
          val(bucket_size), \
          val(slave_tree), \
          file(test_alignment), \
          file(ref_alignment) \
          from toEvaluate

    output:
      set file("*.sp"), \
          file("*.tc"), \
          file("*.col") \
          into scores

    when:
      params.evaluate

     script:
     """
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode sp \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.sp"

       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode tc \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.tc"

       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode column \
            | grep -v "seq1" | grep -v '*' | \
              awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_tree}.col"

    """
}

workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}