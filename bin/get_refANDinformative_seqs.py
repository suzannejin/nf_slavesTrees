#!/usr/bin/env python3

__description__ = '''

Get reference sequences + (N-ref) informative sequences for --> QuanTest2

Warning: You should have T-coffee installed before running this script.

'''


# ==========
# FUNCTIONS 
# ==========

def read_fasta(filename):
    ''' Read a FASTA file.
    
    An ordered ditionary {name : sequence} is returned. 
    '''

    import collections

    d=collections.OrderedDict()

    with open(filename) as f:
        for line in f:
            line=line.strip("\n")
            if line[0]==">":
                name=line[1:]
                d[name]=""
            else:
                d[name]+=line
                
    return(d)


def read_aux(filename):
    ''' Read the lines of a file.

    '''
    l=[]

    with open(filename) as f:
        for line in f:
            line=line.strip("\n")
            l.append(line)

    return(l)


def get_informative_names(fastaname,tree,n):
    ''' Get the name of the informative n sequences from a fasta file <msa or seq>, based on the guide tree <tree> 
    
    An list with the name of the n informative sequences is returned.
    '''

    inf_names=[]
    
    tmp=os.popen("t_coffee -other_pg seq_reformat -in {} -in2 {} -action +regtrim {}".format(fastaname,tree,n)).read()
    tmp=[i for i in tmp.split(">") if i]

    for item in tmp:
        item=[i for i in item.split("\n") if i] # Split into [ name, sequence ]
        name=item[0]
        inf_names.append(name)

    return(inf_names)


def merge_ref_informative_seqs(msa,ref_names,inf_names,n):
    ''' Merge references + (n-ref) informative sequences 

    An ordered dictionary {name : sequence} with n sequences (ref + informative) is returned.
    '''

    import collections

    d=collections.OrderedDict()

    # Reference sequences
    for name in ref_names:
        seq=msa[name]
        d[name]=seq
    
    # Add N-ref informative sequences
    n=len(ref_names)
    for name in inf_names:
        if n>=args.n:
            break
        if name in ref_names:
            continue
        seq=msa[name]
        d[name]=seq
        n+=1
    
    return(d)


# =====
# MAIN
# =====

if __name__ == '__main__':

    import argparse
    import os
    import sys
   

    app = argparse.ArgumentParser(description=__description__)
    app.add_argument('--seq',type=str,help="Sequences in FASTA format. \
                     Note that the informative sequences will be retrieved based on this file (instead of the MSA in order to save time) and the guide tree")
    app.add_argument('--msa',type=str,help="MSA in FASTA format")
    app.add_argument('--names',type=str,help="A file with a list of reference sequence names. \
                     In particular, a total of 3 sequences with known structures are assigned as the references \
                     (whose secondary structure will be predicted using QuanTest2)")
    app.add_argument('--tree',type=str,help="Guide tree to be used to retrieve the informative sequences --> T-coffee +regtrim")
    app.add_argument('--n',type=int,help="Total number of sequences (ref + informative)")
    args = app.parse_args() 

    # Read FASTA
    seqs=read_fasta(args.seq)
    msa=read_fasta(args.msa)

    # Get reference sequence names
    ref_names=read_aux(args.names)
    # Get N informative sequence names
    inf_names=get_informative_names(args.seq,args.tree,args.n)
    # Merge ref + (N-ref) informative sequences
    out=merge_ref_informative_seqs(msa,ref_names,inf_names,args.n)

    # Write output: references + (N-ref) informative sequences
    for name,seq in out.items():
        sys.stdout.write(">"+name+"\n"+seq+"\n")

