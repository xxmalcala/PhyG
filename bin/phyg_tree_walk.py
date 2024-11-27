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

NEEDS SUPPORT FOR UNKNOWN TAXONOMIC AFFILIATIONS!

"""

import glob, sys, time

import numpy as np
import pandas as pd

from collections import defaultdict
from datetime import datetime, timedelta
from ete3 import Tree
from pathlib import Path

# from bin import... MAKE SURE TO ADD THIS AFTER FINISHED TESTING!!!
from bin import phyg_ncbi_taxonomy as tx


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

    # Obnoxiously, this is needed as polytomies are commin in FastTree-based trees...
    t.resolve_polytomy(recursive = True)

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
        leaf_taxonomy: dict,
        taxon_delim: str = '|') -> dict:

    major_clades = ['Bacteria', 'Archaea', 'Opisthokonta','Archaeplastida','Amoebozoa','Discoba','Metamonada','SAR','Orphan']
    clade_sizes = {i:[] for i in ['Prok', 'Opisthokonta','Archaeplastida','Amoebozoa','Discoba','Metamonada','SAR','Orphan']}

    for node in tree.iter_descendants("postorder"):
        mjr_c = []
        node_seqs = [i.name for i in node.get_leaves()]
        node_taxa = list(set([i.partition(taxon_delim)[0] for i in node_seqs]))

        for i in node.get_leaves():
            node_taxa.append(i.name.partition(taxon_delim)[0])
            node_txnmy = leaf_taxonomy[i.name]

            if node_txnmy[0] in ['Bacteria','Archaea']:
                mjr_c.append(node_txnmy[0])

            else:
                mjr_c.append(node_txnmy[1])

        node_taxa = list(set(node_taxa))

        if len(mjr_c) > 1:
            if mjr_c.count('Bacteria') + mjr_c.count('Archaea') == len(mjr_c):
                clade_sizes['Prok'].append((node, len(mjr_c), len(node_taxa), ';'.join(node_seqs)))

            else:
                for clade in major_clades[2:]:
                    if mjr_c.count(clade) == len(mjr_c):
                        taxa = list(set([]))
                        clade_sizes[clade].append((node, len(mjr_c), len(node_taxa), ';'.join(node_seqs)))

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
            largest_clade = max(v, key = lambda x: x[2])

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

    return tree_topo_summary, leaf_taxonomy, tree


def check_many_trees(
        start_time,
        tree_dir: str,
        blen_mode: str = 'median',
        taxon_delim: str = '|',
        query_taxon: str = None,
        verbose: bool = False):

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

    tree_num = 1
    for tree_file in all_tree_files:
        if verbose:
            print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Processing diagnostic phylogeny {tree_num} of {len(all_tree_files)}', end = '\r')
        tfs, lt, t = score_tree(
                        tree_file,
                        reroot_dir,
                        blen_mode,
                        taxon_delim,
                        query_taxon)

        if not tfs:
            continue

        comp_summary += tfs
        tree_num += 1

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Processing diagnostic phylogeny {tree_num} of {len(all_tree_files)}')

    return comp_summary


def paralog_selection(
        tree,
        ortho_dict: dict):

    orthologs_to_keep = []

    for k, v in ortho_dict.items():
        n_dists = []

        max_tax = max([int(i[1]) for i in v])

        tmp = [i for i in v if int(i[1]) == max_tax]

        for i in tmp:
            n_dists.append((i, tree.get_distance(tree.get_leaves_by_name(i[0])[0])))

        top_ortho = min(n_dists, key=lambda x: x[1])[0]

        orthologs_to_keep.append(f'{top_ortho[0]}\t{top_ortho[-1]}')

    return orthologs_to_keep


def ortholog_by_clade(
        tree_file: str,
        reroot_dir: str,
        blen_mode: str,
        taxon_delim: str):

    init_ortho_dict = defaultdict(list)

    cleaner_ortho_dict = {}

    t_topo, leaf_txmy, tree = score_tree(
                                tree_file,
                                reroot_dir,
                                blen_mode,
                                taxon_delim)

    for seq in t_topo:
        d_eval, s_eval, m_eval = 'diff','diff','diff'

        txnmy = leaf_txmy[seq[0]]
        dmn, spr_grp, mjr_cl = seq[3:6]

        if dmn == txnmy[0]:
            d_eval = 'same'

        if spr_grp == txnmy[1]:
            s_eval = 'same'

        if mjr_cl == txnmy[2]:
            m_eval = 'same'

        init_ortho_dict[seq[0]].append([seq[1],seq[7], seq[-2], d_eval, s_eval, m_eval, seq[8]])

    for k, v in init_ortho_dict.items():
        seqs_to_keep = [i for i in v if i.count('same') == 3]

        if not seqs_to_keep:
            seqs_to_keep = [i for i in v if i.count('same') == 2]

        if not seqs_to_keep:
            seqs_to_keep = [i for i in v if i.count('same') == 1]

        else:
            seqs_to_keep = [i for i in v]

        cleaner_ortho_dict[k] = seqs_to_keep

    return list(init_ortho_dict.values()), paralog_selection(tree, cleaner_ortho_dict)


def eval_many_clades(
    project_name: str,
    tree_dir: str,
    blen_mode: str = 'median',
    taxon_delim: str = '|'):

    orthologs_out = f'{project_name}.Ortholog_Selections.txt'

    reroot_dir = f'{tree_dir.rstrip("/")}_ReRooted/'
    out_clade_dir = f'{reroot_dir.rpartition("_")[0]}_Clade_Size_Summary/'

    Path(reroot_dir).mkdir(parents = True, exist_ok = True)
    Path(out_clade_dir).mkdir(parents = True, exist_ok = True)

    ortholog_selections = []

    all_tree_files = glob.glob(f'{tree_dir}/*.tre*')

    if not all_tree_files:
        all_tree_files = glob.glob(f'{tree_dir}/*.nwk*')

    if not all_tree_files:
        all_tree_files = glob.glob(f'{tree_dir}/*newick')

    if not all_tree_files:
        print('ERROR: Unable to find phylogenetic trees ending with: ".tre*", ".nwk", nor ".newick"\n')
        print('       Please adjust the ending extension of your phylogenetic trees.')
        sys.exit(1)


        """
        For now just goes by major clade, could be more specific (e.g., minor clade)...
        """

    for tree_file in all_tree_files:

        out_tsv = f'{out_clade_dir}{tree_file.rpartition("/")[-1].partition(".")[0]}.Clade_Sizes.tsv'

        topo_summary, ologs = ortholog_by_clade(
                                tree_file,
                                reroot_dir,
                                blen_mode,
                                taxon_delim)

        topology_info = ['\t'.join(i) for j in topo_summary for i in j]
        ortholog_selections += ologs


        with open(out_tsv, 'w+') as w:
            w.write('Sequence\tNumber_Unique_Sister_Taxa\tBranch_Length\tSame_Domain\t' \
                'Same_Supergroup\tSame_Major_Clade\tSister_Taxa')
            w.write('\n'.join(topology_info))

    with open(orthologs_out, 'w+') as w:
        w.write('Sequence\tSister_Taxa\n')
        w.write('\n'.join(ortholog_selections))


if __name__ == '__main__':
    try:
        project_name = sys.argv[1]
        tree_dir = sys.argv[2]
    except:
        print('Usage:\n\n    python3 phyg_tree_walk.py [PROJECT-NAME] [TREE-DIRECTORY]\n')
        sys.exit()

    eval_many_clades(
        project_name,
        tree_dir)
