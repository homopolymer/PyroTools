#!/usr/bin/env python
"""
SYNOPSIS

    TODO plot_taxonomy_tree [-h,--help] [-v,--verbose] <NEWICK> <TAXONOMY>

DESCRIPTION

    TODO This provides a tree view of OTU taxonomy.

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

import os
import tempfile
from contextlib import contextmanager
from cStringIO import StringIO
from Bio import Phylo
from ete2 import Phyloxml,faces,AttrFace,TreeStyle
from sets import Set

# ref: http://stackoverflow.com/questions/11892623/python-stringio-and-compatibility-with-with-statement-context-manager
@contextmanager
def tempinput(data):
    temp = tempfile.NamedTemporaryFile(delete=False)
    temp.write(data)
    temp.close()
    yield temp.name
    os.unlink(temp.name)


def tree_layout(node):
    if not node.is_leaf():
        #faces.add_face_to_node(AttrFace("mylabel"), node, column=1, position="branch-top")
        faces.add_face_to_node(AttrFace("mylabel"), node, column=1)
    else:
        if hasattr(node,'color'):
            faces.add_face_to_node(AttrFace("name",fgcolor=node.color), node, column=1)
        else:
            faces.add_face_to_node(AttrFace("name"), node, column=1)

def get_taxonomy(fn):
    tax = {}
    with open(fn) as f:
        for line in f:
            item = line.split()
            tax[item[0]] = {'kingdom':item[1][:-1], \
                            'phylum': item[2][:-1], \
                            'class':  item[3][:-1], \
                            'order':  item[4][:-1], \
                            'family': item[5][:-1], \
                            'genus':  item[6][:-1], \
                            'species':item[7][:-1]}
    return tax


def add_taxonomy_for_internal_branch(tree,tax):
    iid = 0
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.add_feature("id",node.name)
        else:
            node.add_feature("id","i"+str(iid))
            tax["i"+str(iid)]={}
            for t in ['kingdom','phylum','class','order','family','genus','species']:
                ts = Set([tax[x.id][t] for x in node.get_children()])
                if len(ts)==1:
                    tax[node.id][t]=ts.pop()
                else:
                    tax[node.id][t]="%s__"%t[0]
            iid = iid + 1
    return tree



def add_node_label(tree,tax):
    for node in tree.traverse():
        if not node.is_leaf():
            label = "null"
            for t in ['kingdom','phylum','class','order','family','genus','species']:
                if len(tax[node.id][t])>3:
                    label = tax[node.id][t]
            node.add_feature("mylabel",label)
    return tree

   
def main ():

    global options, args

    if options.verbose: print time.asctime(),
    if options.verbose: print "load and parse newick file"
    # TODO: read newick file
    tree = Phylo.read(args[0],'newick') 
    # TODO: convert newick to phyloxml
    treeXML = StringIO()
    Phylo.write(tree,treeXML,'phyloxml')
    # TODO: read phyloxml as ete object
    hPhylotree = Phyloxml()
    with tempinput(treeXML.getvalue()) as tempfilename:   
        hPhylotree.build_from_file(tempfilename)
    # TODO: get the tree
    tree2 = hPhylotree.get_phylogeny()[0]
       
    if options.verbose: print time.asctime(),
    if options.verbose: print "load and parse taxonomy file"
    # TODO: read taxonomy file
    tax = get_taxonomy(args[1])
    # TODO: refine taxonomy annotation of internal node
    tree2 = add_taxonomy_for_internal_branch(tree2,tax)
    # TODO: refine tree node label
    #tree2 = add_node_label(tree2,tax)
    for node in tree2.traverse():
        if not node.is_leaf():
            label = "null"
            for t in ['kingdom','phylum','class','order','family','genus','species']:
                if len(tax[node.id][t])>3:
                    label = tax[node.id][t]
            node.add_feature("mylabel",label)

    # TODO: add node depth
    depth={}
    if options.depth:
        with open(options.depth) as f:
            for line in f:
                (id,dep) = line.split()
                depth[id] = float(dep)
   
    # TODO: add color attribute
    if options.depth:
        for node in tree2.iter_leaves():
            if depth[node.id] >= 10 and depth[node.id] < 100:
                node.add_feature("color","#D8BFD8")
            elif depth[node.id] >= 100 and depth[node.id] < 1000:
                node.add_feature("color","#DDA0DD")
            elif depth[node.id] >= 1000 and depth[node.id] < 5000:
                node.add_feature("color","#EE82EE")
            elif depth[node.id] >= 5000:
                node.add_feature("color","#DA70D6")
            else:
                node.add_feature("color","#E6E6FA")
        
    # TODO: set tree style
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = tree_layout
    # TODO: show tree2
    #tree2.show(tree_style=ts)
    tree2.render(args[0]+".png",dpi=2048,tree_style=ts)
    

if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'], version='$Id$')
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option ('-d', '--depth', help='sequencing depth', type='str')
        (options, args) = parser.parse_args()
        if len(args) < 2:
            parser.error ('missing argument')
        main()
        if options.verbose: print 'TOTAL TIME IN MINUTES:',
        if options.verbose: print (time.time() - start_time) / 60.0
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

