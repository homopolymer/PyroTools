#!/usr/bin/env python
"""
SYNOPSIS

    TODO pick_representative_otus.py [-h,--help] [-v,--verbose] [-c,--clade <LEVEL>] <TAXONOMY> <DEPTH> <FASTA>

DESCRIPTION

    TODO This script picks out the representative OTU(s) 
    within a specified clade level.

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

AUTHOR

    TODO: Feng Zeng <zengfeng@xmu.edu.cn>

LICENSE

    This script is in the public domain, free from copyrights or restrictions.

VERSION

    $Id$
"""

import sys, os, traceback, optparse
import time
import re

import operator
from ete2 import Phyloxml,faces,AttrFace,TreeStyle,Tree
from sets import Set
import numpy as np
from skbio.alignment import local_pairwise_align_ssw, Alignment


def ssw_similarity(seq_a,seq_b):
    a = local_pairwise_align_ssw(seq_a,seq_b)
    return 1-a.distances()[0][1]


def tree_layout(node):
    if not node.is_leaf() and node.name:
        #faces.add_face_to_node(AttrFace("mylabel"), node, column=1, position="branch-top")
        faces.add_face_to_node(AttrFace("name"), node, column=1)
    if node.is_leaf():
        faces.add_face_to_node(AttrFace("depth",fgcolor="red"), node, column=1)


def get_gg_taxonomy(fn):
    tax = {}
    with open(fn) as f:
        for line in f:
            (id,annot) = line.split("\t")
            item = annot.split(";")
            for i in range(len(item)): item[i] = item[i].strip()
            tax[id] = {'kingdom':item[0], \
                       'phylum': item[1], \
                       'class':  item[2], \
                       'order':  item[3], \
                       'family': item[4], \
                       'genus':  item[5], \
                       'species':item[6]}
    return tax


def main ():

    global options, args

    # TODO: read and parse taxonomy file
    if options.verbose: print >> sys.stderr, "[",time.asctime(),"]",
    if options.verbose: print >> sys.stderr, "parse taxonomy file"
    tax = get_gg_taxonomy(args[0])
    
    # TODO: construct the taxonomy tree
    if options.verbose: print >> sys.stderr, "[",time.asctime(),"]",
    if options.verbose: print >> sys.stderr, "build taxonomy tree"
    tree = Tree()
    tree_node = {}
    for x in tax:
        n = tree.search_nodes(name=tax[x]['kingdom'])
        if len(n)==0:
            n = tree.add_child(name=tax[x]['kingdom'])
            tree_node[tax[x]['kingdom']] = n

    for x in tax:
        y_name = tax[x]['kingdom']
        for level in ['phylum','class','order','family','genus','species']:
            if len(tax[x][level])>3:
                x_name = y_name+"_"+tax[x][level]
                if x_name not in tree_node:
                    n = tree_node[y_name].add_child(name=tax[x][level])
                    tree_node[x_name] = n
                y_name = x_name
            else:
                n = tree_node[y_name].add_child(name=x)
                tree_node[x] = n
                break
            if level=='species':
                n = tree_node[y_name].add_child(name=x)
                tree_node[x] = n

    # TODO: read and parse depth file
    if options.verbose: print >> sys.stderr, "[",time.asctime(),"]",
    if options.verbose: print >> sys.stderr, "parse depth file"
    depth = {}
    with open(args[1]) as f:
        for line in f:
            item = line.split()
            depth[item[0]] = float(item[1])

    # TODO: read and parse fasta file
    from skbio.parse.sequences import parse_fasta
    if options.verbose: print >> sys.stderr, "[",time.asctime(),"]",
    if options.verbose: print >> sys.stderr, "parse fasta file"
    fasta = {}
    for id,seq in parse_fasta(args[2]):
        fasta[id] = seq

    # TODO: prune the tree to remove empty leaves 
    retain_nodes = Set()
    for node in tree.iter_leaves():
        if depth[node.name]>0 and len(tax[node.name][options.clade])>3:
            retain_nodes.add(node)
            retain_nodes.update(Set(node.get_ancestors()))

    tree.prune(retain_nodes)

    # TODO: collect clade statistics
    clade_depth = {}
    for node in tree.traverse():
        if node.name.startswith(options.clade[0]):
            d = 0
            for n in node.get_leaves():
                if d < depth[n.name]:
                    d = depth[n.name]
            clade_depth[node.name] = d
    
    # TODO: get the threshold
    X = np.array(sorted(depth.values(),reverse=True))
    rX = X.cumsum()*100/X.sum()
    T = np.interp(options.Npercent,rX,X)

    # TODO: filter out clades below the threshold
    clade_retain = Set()
    for c,d in clade_depth.items():
        if d>=T:
            clade_retain.add(c) 

    # TODO: prune tree again
    retain_nodes = Set()
    for node in tree.get_leaves():
        if tax[node.name][options.clade] in clade_retain:
            retain_nodes.add(node)
            retain_nodes.update(Set(node.get_ancestors()))

    tree.prune(retain_nodes)

    # TODO: pick out representative OTUs for a clade
    retain_nodes = Set()
    for taxon in clade_retain:
        t = tree & taxon
        L = {n:depth[n.name] for n in t.get_leaves()}
        l = sorted(L.items(),key=operator.itemgetter(1),reverse=True)
        taxon_otus = Set()
        for n, val in l:
            if val<T: continue
            if not taxon_otus:
                taxon_otus.add(n)
            else:
                max_sim = 0
                for y in taxon_otus:
                    sim = ssw_similarity(fasta[n.name],fasta[y.name])
                    if max_sim < sim: max_sim = sim
                if max_sim < options.Tsim:
                    taxon_otus.add(n)
        retain_nodes.update(taxon_otus)
        for n in taxon_otus:
            retain_nodes.update(Set(n.get_ancestors()))
        
    # TODO: pick out a representative OTUs for a clade
    #retain_nodes = Set()
    #for taxon in clade_retain:
    #    t = tree & taxon
    #    L = {n:depth[n.name] for n in t.get_leaves()}
    #    l = sorted(L.items(),key=operator.itemgetter(1),reverse=True)
    #    retain_nodes.add(l[0][0])
    #    retain_nodes.update(Set(l[0][0].get_ancestors()))
    #
    tree.prune(retain_nodes)
 
    # TODO: output representative OTUs
    for node in tree.get_leaves():
        print node.name

    #
    for node in tree.get_leaves():
        node.add_feature("depth",depth[node.name])

    if options.png:
        ts = TreeStyle()
        ts.layout_fn = tree_layout
        tree.render(options.png,dpi=1024,tree_style=ts)
                
        
    

if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'], version='$Id$')
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option ('-c', '--clade', default='phylum', help='specify clade level', type='str')
        parser.add_option ('-n', '--Npercent', default=95, help='specify depth percentage level for picking OTUs, like the metric N50', type='float')
        parser.add_option ('-s', '--Tsim', default=0.85, help='specify sequence similarity threshold', type='float')
        parser.add_option ('-p', '--png', help='draw the tree for representative OTUs, specify filename', type='str')
        (options, args) = parser.parse_args()
        if len(args) < 3:
            parser.error ('missing argument')
        main()
        if options.verbose: print >> sys.stderr, 'TOTAL TIME IN MINUTES:',
        if options.verbose: print >> sys.stderr, (time.time() - start_time) / 60.0
        sys.exit(0)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)

