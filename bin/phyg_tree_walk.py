#!/usr/bin/env python3

"""
Phylogenetic tree-walking based approach for assessing sister relationships
in single-gene phylogenetic trees.

Will be VERY simple and include:
- Largest clades, based on katz-lab taxonomy concepts (Major-Minor Clades) OR
  user-defined table of clades.
- Single-sister relationships
- use NCBI's taxonomy if no clade information given
- - this will need name delimiters!
- - can try naively with the initial split(from the delimiter) and if that fails
- - consistently, then default to the next position?

Outputs, big table, summary table based on single-sister relationships and short
branches (either median or mean of all branches)...

"""

import glob, sys

import numpy as np
import pandas as pd

from collections import defaultdict
from ete3 import Tree
from pathlib import Path

# from bin import... MAKE SURE TO ADD THIS AFTER FINISHED TESTING!!!
import phyg_ncbi_taxonomy as tx


"""
does not support custom clade assignments yet...
"""


def parse_tree_file(
        tree_file: str,
        taxon_delim: str = '|',
        custom_clade = False,
        clade_info: str = None,
        clade_level: str = None) -> list:

    "Clade info will be the custom table of clades..."

    t = Tree(tree_file)

    leaf_taxonomy = {}
    leaf_names = [node.name for node in t.get_leaves()]

    for leaf in leaf_names:
        taxon_name = leaf.partition(taxon_delim)[0]

        txnmy = tx.search_ncbi_lineage(leaf, False)

        leaf_taxonomy[leaf] = txnmy
        leaf_taxonomy[taxon_name] = txnmy

    if len(set([v[1] for v in leaf_taxonomy.values()])) < 2:
        print(f'Ignoring {tree_file} as it contains fewer than 2 "major clades"')
        return None

    else:
        return t, leaf_taxonomy


def save_tree(
        tree,
        tree_file: str,
        reroot_dir: str) -> None:

    out_tree = f'{reroot_dir}{tree_file.rpartition("/")[-1].rpartition(".")[0]}.ReRooted.nwk'

    tree.write(format = 1, outfile = out_tree)


def eval_clade_sizes(
        tree,
        leaf_taxonomy: dict) -> dict:

    major_clades = ['Bacteria', 'Archaea', 'Opisthokonta','Archaeplastida','Amoebozoa','Discoba','Metamonada','SAR','Orphan']
    clade_sizes = {i:[] for i in ['Prok', 'Opisthokonta','Archaeplastida','Amoebozoa','Discoba','Metamonada','SAR','Orphan']}

    for node in tree.iter_descendants("postorder"):
        mjr_c = []

        for i in node.get_leaves():
            node_txnmy = leaf_taxonomy[i.name]

            if node_txnmy[0] in ['Bacteria','Archaea']:
                mjr_c.append(node_txnmy[0])

            else:
                mjr_c.append(node_txnmy[1])

        if len(mjr_c) > 1:
            if mjr_c.count('Bacteria') + mjr_c.count('Archaea') >= (len(mjr_c) - 1):
                clade_sizes['Prok'].append((node, len(mjr_c)))

            else:
                for clade in major_clades[2:]:
                    if mjr_c.count(clade) >= (len(mjr_c) - 1):
                        clade_sizes[clade].append((node, len(mjr_c)))

    return clade_sizes


def adjust_tree_root(
        tree,
        leaf_taxonomy: dict):
    """
    need to support custom clades!!!!
    """

    tree.set_outgroup(tree.get_midpoint_outgroup())
    clade_sizes = eval_clade_sizes(tree, leaf_taxonomy)

    for k, v in clade_sizes.items():
        if v:
            largest_clade = max(v, key = lambda x: x[-1])
            if largest_clade[0]:
                tree.set_outgroup(largest_clade[0])
                break
        else:
            continue

    return tree


def check_same_taxon(
        node,
        taxon_delim: str = '|',
        reps: int = 0):

    taxon_node_names = list(set([leaf.name.partition(taxon_delim)[0] for leaf in node.up.get_leaves()]))

    if len(taxon_node_names) == 1:
        parent = node.up
        return parent, True, reps + 1

    else:
        return node.up, False, reps


def eval_sister_relationships(
        node,
        sister_seqs: list,
        leaf_taxonomy: dict,
        taxon_delim: str,
        thresh_blen: float):

    sister_taxa = list(set([i.partition(taxon_delim)[0] for i in sister_seqs]))

    dmn_c = []
    spr_c = []
    mjr_c = []
    mnr_c = []

    for i in sister_taxa:
        dmn_c.append(leaf_taxonomy[i][0])
        spr_c.append(leaf_taxonomy[i][1])
        mjr_c.append(leaf_taxonomy[i][2])
        mnr_c.append(leaf_taxonomy[i][3])

    if len(set(dmn_c)) > 1:
        dmn_eval = 'No-clear-domain'
    else:
        dmn_eval = dmn_c[0]

    if len(set(spr_c)) > 1:
        spr_eval = 'No-clear-major-group'
    else:
        spr_eval = spr_c[0]

    if len(set(mjr_c)) > 1:
        mjr_eval = 'No-clear-major-clade'
    else:
        mjr_eval = mjr_c[0]

    if len(set(mnr_c)) > 1:
        mnr_eval = 'No-clear-minor-clade'
    else:
        mnr_eval = mnr_c[0]

    if node.dist <= thresh_blen:
        qblen_type = 'short'
    else:
        qblen_type = 'long'

    if len(sister_taxa) <= 20:
        sis_tax_list = ';'.join(sister_taxa)
        sis_seq_list = ';'.join(sister_seqs)
    else:
        sis_tax_list = 'More-than-20'
        sis_seq_list = 'Too-many-seqs'

    node_summary = [
        dmn_eval,
        spr_eval,
        mjr_eval,
        mnr_eval,
        f'{len(set(sister_taxa))}',
        sis_tax_list,
        sis_seq_list,
        f'{node.dist:.4f}',
        f'{thresh_blen:.4f}',
        qblen_type]

    return node_summary


def check_sisters(
        tree,
        leaf_taxonomy,
        blen_mode,
        taxon_delim: str = '|',
        query_taxon: str = None):

    tree_phylo_sister_summary = []

    tip_ancest_dists = [i.dist for i in tree.get_leaves()]

    if blen_mode == 'average':
        thresh_blen = np.mean(tip_ancest_dists)

    elif blen_mode == 'median':
        thresh_blen = np.median(tip_ancest_dists)

    for node in tree.get_leaves():
        if query_taxon and query_taxon not in node.name:
            continue

        leaf_query_taxon_eval = [node, True, 0]
        leaf_query_taxon = node.name.partition(taxon_delim)[0]

        while leaf_query_taxon_eval[1]:
            leaf_query_taxon_eval = check_same_taxon(
                                        leaf_query_taxon_eval[0],
                                        taxon_delim,
                                        leaf_query_taxon_eval[2])

        node_seqs = [taxon.name for taxon in leaf_query_taxon_eval[0].get_leaves()]

        sister_seqs = [i for i in node_seqs if i.partition(taxon_delim)[0] != leaf_query_taxon]

        tmp_summary = eval_sister_relationships(
                        node,
                        sister_seqs,
                        leaf_taxonomy,
                        taxon_delim,
                        thresh_blen)

        node_summary = [
            leaf_query_taxon,
            node.name,
            node.name.rpartition(taxon_delim)[-1]
            ] + tmp_summary

        tree_phylo_sister_summary.append(node_summary)

    return tree_phylo_sister_summary


def score_tree(
    tree_file: str,
    reroot_dir: str,
    blen_mode: str,
    taxon_delim: str,
    query_taxon: str = None) -> list:

    tree, leaf_taxonomy = parse_tree_file(tree_file)

    if not tree:
        return None

    tree = adjust_tree_root(tree, leaf_taxonomy)

    save_tree(tree, tree_file, reroot_dir)

    tree_topo_summary = check_sisters(
                            tree,
                            leaf_taxonomy,
                            blen_mode,
                            taxon_delim,
                            query_taxon)

    return tree_topo_summary


def check_many_trees(
        tree_dir: str,
        blen_mode: str = 'median',
        taxon_delim: str = '|',
        query_taxon: str = None):

    reroot_dir = f'{tree_dir.rstrip("/")}_ReRooted/'

    Path(reroot_dir).mkdir(parents = True, exist_ok = True)

    comp_summary = []
    all_tree_files = glob.glob(f'{tree_dir}/*.tre*')

    if not all_tree_files:
        all_tree_files = glob.glob(f'{tree_dir}/*.nwk*')

    if not all_tree_files:
        all_tree_files = glob.glob(f'{tree_dir}/*newick')

    if not all_tree_files:
        print('ERROR: Unable to find phylogenetic trees ending with: ".tre*", ".nwk", nor ".newick"\n')
        print('       Please adjust the ending extension of your phylogenetic trees.')
        sys.exit(1)

    for tree_file in all_tree_files:
        print(tree_file)
        tree_file_summary = score_tree(
                                tree_file,
                                reroot_dir,
                                blen_mode,
                                taxon_delim,
                                query_taxon)

        if not tree_file_summary:
            continue

        comp_summary += tree_file_summary

    return comp_summary
