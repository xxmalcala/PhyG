#!/usr/bin/env python3

import glob, shutil, subprocess, sys, time

from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

from Bio import SeqIO

import pandas as pd

# from bin import... MAKE SURE TO ADD THIS AFTER FINISHED TESTING!!!
from bin import phyg_add_taxa as at
from bin import phyg_ncbi_taxonomy as tx
from bin import phyg_tree_walk as pt


"""
Notes for self:

pull prep_peptide_ngrams and ngram_from_seq from phyg_msa_prep?

need taxonomy dictionary for the contamination database or update the taxon names
in the database and use NCBITaxa from ete3 to do the heavy lifting (might be more work)

Just does the basics for now... need to incorporate the 'one-offs' for the TIdeS
decontamination tests.
"""

# MOVE THIS OR UPDATE SLIGHTLY?
def prep_query_data(
        start_time,
        fasta_file: str,
        query_taxon: str,
        contam_dir: str,
        diag_db: str,
        gen_code: str = '1',
        taxon_code: str = None,
        transcripts: bool = True,
        genbank: bool = False,
        prokaryotic: bool = False,
        threads: int = 4,
        verbose: bool = False) -> None:

    if glob.glob(f'{contam_dir}*AA.fasta'):
        # print('query data is already prepared')
        out_pep_fas = glob.glob(f'{contam_dir}*AA.fasta')[0]

        return contam_dir, out_pep_fas

    elif transcripts:
        at.prep_transcriptomes(
                start_time,
                fasta_file,
                contam_dir,
                query_taxon,
                taxon_code,
                diag_db,
                min_hit_cover = 10,
                gen_code = gen_code,
                threads = threads,
                verbose = verbose)
    else:
        at.prep_wgs(
                start_time,
                fasta_file,
                contam_dir,
                query_taxon,
                taxon_code,
                diag_db,
                min_hit_cover = 10,
                gen_code = gen_code,
                genbank = genbank,
                prokaryotic = prokaryotic,
                threads = threads,
                verbose = verbose)

    out_pep_fas = glob.glob(f'{contam_dir}*AA.fasta')[0]

    return contam_dir, out_pep_fas


def prep_query_peptides(
        out_dir: str,
        query_orf_fasta: str) -> list:

    prepped_orf_dir = f'{out_dir}Phylo_Based_Contam/Query_ORFs/'
    Path(prepped_orf_dir).mkdir(parents = True, exist_ok = True)

    og_db = defaultdict(list)

    phylo_file_dict = defaultdict(list)

    for i in SeqIO.parse(query_orf_fasta,'fasta'):
        og_db[i.id.rpartition("|")[-1]].append(i)

    for k, v in og_db.items():
        SeqIO.write(v, f'{prepped_orf_dir}{k}.Query_Seqs.AA.fasta','fasta')
        phylo_file_dict[k].append(f'{prepped_orf_dir}{k}.Query_Seqs.AA.fasta')

    return prepped_orf_dir, phylo_file_dict


def mafft_align(
        query_msa_dir: str,
        query_taxon: str,
        diag_og: str,
        query_seq_fasta: str,
        ref_msa_fasta: str,
        threads: int = 4) -> str:

    out_query_msa = f'{query_msa_dir}{diag_og}.Query_MSA.AA.fasta'

    mafft_cmd = f'mafft --quiet ' \
                f'--thread {threads} ' \
                '--auto ' \
                f'--add {query_seq_fasta} ' \
                f'--keeplength {ref_msa_fasta} ' \
                f'> {out_query_msa}'

    mafft_rslt = subprocess.Popen(
                    mafft_cmd,
                    stdout = subprocess.PIPE,
                    stderr = subprocess.PIPE,
                    universal_newlines = True,
                    shell = True)

    stdout, sderr = mafft_rslt.communicate()

    tmp = [i for i in SeqIO.parse(out_query_msa,'fasta') if query_taxon in i.id]
    SeqIO.write(tmp, out_query_msa, 'fasta')

    return out_query_msa


def convert_jplace_newick(
        outdir: str,
        jplace_file: str) -> None:

    gappa_cmd = f'gappa examine graft ' \
                '--fully-resolve ' \
                f'--jplace-path {jplace_file} ' \
                f'--out-dir {outdir}'

    gappa_rslt = subprocess.run(
                    gappa_cmd.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)


def run_epa_ng(
        updated_tree_dir: str,
        ref_msa: str,
        ref_tree: str,
        query_msa: str,
        threads = 4):

    out_jplace = f'{updated_tree_dir}{query_msa.rpartition("/")[-1].partition(".")[0]}.Updated.jplace'

    epa_ng_cmd = f'epa-ng ' \
                    f'-T {threads} ' \
                    f'--ref-msa {ref_msa} ' \
                    f'--tree {ref_tree} ' \
                    f'--query {query_msa} ' \
                    f'-w {updated_tree_dir} ' \
                    '--model LG+I+G4'

    # print(epa_ng_cmd)

    epa_ng_rslt = subprocess.run(
                    epa_ng_cmd.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)

    shutil.copy2(f'{updated_tree_dir}epa_result.jplace',f'{out_jplace}')

    convert_jplace_newick(updated_tree_dir, out_jplace)

    Path.unlink(Path.cwd() / f'{updated_tree_dir}epa_info.log')
    Path.unlink(Path.cwd() / f'{updated_tree_dir}epa_result.jplace')
    Path.unlink(Path.cwd() / out_jplace)


def prep_peptide_kmers(
        gf_seqs: list,
        kmer_len: int = 5,
        overlap: bool = True) -> dict:

    peptide_kmers = {i.id:kmer_from_seq(f'{i.seq}', kmer_len, overlap = overlap) for i in gf_seqs}

    return peptide_kmers


def kmer_from_seq(
        seq: str,
        kmer_len: int,
        overlap: bool) -> str:

    if overlap:
        return ' '.join([seq[n:n+kmer_len] for n in range(len(seq) - kmer_len + 1)])

    else:
        return ' '.join([seq[n:n+ngram_len] for n in range(0, len(seq) - kmer_len + 1, kmer_len)])


def knn_cluster(
        start_time,
        query_orf_fasta: str,
        og_reference_fasta: str,
        query_taxon: str,
        kmer_len: int = 5,
        overlap: bool = True,
        clust_thresh: int = 3,
        eval_thresh: float = 0.6,
        verbose: bool = False):

    from sklearn.feature_extraction.text import TfidfVectorizer
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.neighbors import kneighbors_graph

    out_data = []
    gf_num = 0

    og_dict = defaultdict(list)
    contam_taxonomy = {query_taxon: None}

    # Prep the term frequency-inverse document frequency vectorizer
    tfidf = TfidfVectorizer(
                ngram_range = (1,5),
                min_df = 5)

    # Just keep gene families with the query data in them
    for i in SeqIO.parse(query_orf_fasta,'fasta'):
        og_dict[i.id.rpartition("|")[-1]].append(i)

    for i in SeqIO.parse(og_reference_fasta,'fasta'):
        if i.id.rpartition("|")[-1] in og_dict.keys():
            taxon_name = i.id.partition("|")[0]
            og_dict[i.id.rpartition("|")[-1]].append(i)
            contam_taxonomy[taxon_name] = None

    # Capture the reduced taxonomy for all taxa
    for k in contam_taxonomy.keys():
        contam_taxonomy[k] = tx.search_ncbi_lineage(k, False)

    # KNN-graphing, hierarchical clustering, and summary for each gene family with the query taxon
    for k, v in og_dict.items():
        gf_num += 1
        if verbose:
            print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Processing diagnostic gene family {gf_num} of {len(og_dict)}', end = '\r')
        query_cluster_names = []
        eval_clusters = defaultdict(list)
        query_clusters = defaultdict(dict)

        # Generate list of peptide kmers and associated sequence labels
        peptide_kmers = prep_peptide_kmers(v, kmer_len = kmer_len, overlap = overlap)
        peptide_labels = list(peptide_kmers.keys())

        # Generate the TF-IDF array
        peptide_tfidf = tfidf.fit_transform(peptide_kmers.values()).toarray()

        # Build KNN-graph for the gene family's "vocabulary distances"
        knn_graph = kneighbors_graph(
                        peptide_tfidf,
                        10,
                        include_self = False)

        # Make clusters of peptides based on the KNN-graph of "vocab distances"
        agg_mdl = AgglomerativeClustering(
                    distance_threshold = clust_thresh,
                    n_clusters = None)

        agg_mdl = agg_mdl.fit(knn_graph.toarray())

        # Only really care about clusters within the gene family that possess the
        # query taxon
        for n in range(len(agg_mdl.labels_)):
            s = peptide_labels[n]
            c = agg_mdl.labels_[n]
            eval_clusters[c].append(s)

            if query_taxon in s:
                query_cluster_names.append(c)

        # Parse all clusters with query taxa in them
        for n in list(set(query_cluster_names)):
            s = eval_clusters[n]

            # For Katz-lab naming scheme! Needs to be updated!
            # Can use a simple dictionary if using a given database (save some time!)
            q_seqs = [i for i in s if query_taxon in i]
            c_seqs = list(set([i for i in s if query_taxon not in i]))

            # Check to see if non-query taxa are in the cluster
            if c_seqs:
                query_cluster_summary = eval_query_clusters(
                                            contam_taxonomy,
                                            q_seqs,
                                            c_seqs,
                                            eval_thresh)

                out_data += query_cluster_summary



    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Processing diagnostic gene family {gf_num-1} of {len(og_dict)}')

    return out_data


def eval_query_clusters(
        contam_taxonomy: dict,
        query_seqs: list,
        cluster_seqs: list,
        eval_thresh: float = 0.6) -> list:

    summary_data = []

    c_seq_txnmy = [list(contam_taxonomy[i.partition("|")[0]])+[i] for i in cluster_seqs]

    for q in query_seqs:
        domain, spr_grp, mjr, mnr = [], [], [], []

        for i in c_seq_txnmy:
            domain.append(i[0])
            spr_grp.append(i[1])
            mjr.append(i[2])
            mnr.append(i[3])

        # Only care about those clades that are well represented in the cluster
        # 30% seems like a good arbitrary threshold
        domain_eval = [i for i in list(set(domain)) if domain.count(i) > len(domain)*eval_thresh]
        spr_grp_eval = [i for i in list(set(spr_grp)) if spr_grp.count(i) > len(spr_grp)*eval_thresh]
        mjr_eval = [i for i in list(set(mjr)) if mjr.count(i) > len(mjr)*eval_thresh]
        mnr_eval = [i for i in list(set(mnr)) if mnr.count(i) > len(mjr)*eval_thresh]

        # Skip clusters with no clear super-group ties (probably "junk" clusters)
        if not spr_grp_eval:
            continue

        else:
            for grp in spr_grp_eval:
                grp_contams = [i for i in c_seq_txnmy if i[1] == grp]

                for i in grp_contams:
                    x = i
                    if not mjr_eval:
                        x = x[:2]+['No-clear-major-clade']+x[3:]

                    if not mnr_eval:
                        x = x[:3]+['No-clear-minor-clade']+x[4:]

                    x = '\t'.join(x)
                    summary_data.append(f'{q}\t{x}\t{i[-1].partition("|")[0]}')
    if summary_data:
        return summary_data

    else:
        return []


def phylo_contam(
        start_time,
        contam_dir: str,
        query_fas: str,
        ref_msa: str,
        ref_trees: str,
        query_taxon: str,
        blen_mode: str = 'median',
        threads: int = 4,
        verbose: bool = False):

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Starting Phylogeny-based By-Catch Check')

    phylo_dir, phylo_file_dict = prep_query_peptides(
                                    contam_dir,
                                    query_fas)

    query_msa_dir = phylo_dir.replace("_ORFs","_MSAs")
    updated_tree_dir = f'{contam_dir}Phylo_Based_Contam/Updated_Trees/'

    Path(query_msa_dir).mkdir(parents = True, exist_ok = True)
    Path(updated_tree_dir).mkdir(parents = True, exist_ok = True)

    for msa in glob.glob(f'{ref_msa}/*fasta'):
        og = msa.rpartition("/")[-1].partition(".")[0]

        if og in phylo_file_dict:
            phylo_file_dict[og].append(msa)

    for tree in glob.glob(f'{ref_trees}/*nwk'):
        og = tree.rpartition("/")[-1].partition(".")[0]

        if og in phylo_file_dict:
            phylo_file_dict[og].append(tree)

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Adding Query Data to Reference MSAs')

    for k, v in phylo_file_dict.items():
        query_msa = mafft_align(
                        query_msa_dir,
                        query_taxon,
                        k,
                        v[0],
                        v[1],
                        threads = threads)

        phylo_file_dict[k].append(query_msa)

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Adding Query Data to Reference Trees')

    for k, v in phylo_file_dict.items():
        if len(v) == 4:
            run_epa_ng(
                updated_tree_dir,
                v[1],
                v[2],
                v[3],
                threads = threads)

    contam_trees_summary = pt.check_many_trees(
                                start_time,
                                updated_tree_dir,
                                blen_mode,
                                query_taxon = query_taxon,
                                verbose = verbose)

    phylo_contam_eval_trees(
        start_time,
        contam_trees_summary,
        query_taxon,
        updated_tree_dir.rpartition("Updated_Trees/")[0],
        contam_dir,
        blen_mode,
        verbose)


def knn_contam_eval(
        start_time,
        out_dir: str,
        query_orf_fasta: str,
        og_reference_fasta: str,
        query_taxon: str,
        kmer_len: int = 5,
        overlap: bool = True,
        clust_thresh: int = 3,
        eval_thresh: float = 0.6,
        verbose: bool = False):

    print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Starting KNN-based By-Catch Check')

    summary_dir = f'{out_dir}AlignFree_Contam/'
    Path(summary_dir).mkdir(parents = True, exist_ok = True)

    out_base_tsv = f'{summary_dir}{query_taxon}.PhyG_AlignFree_Contam.'

    out_tsv = f'{out_base_tsv}tsv'

    knn_summary_data = knn_cluster(
                    start_time,
                    query_orf_fasta,
                    og_reference_fasta,
                    query_taxon,
                    kmer_len,
                    overlap,
                    clust_thresh,
                    eval_thresh,
                    verbose)

    with open(out_tsv,'w+') as w:
        w.write('Query_Seq\tDomain\tMajor_Group\tMajor_Clade\tMinor_Clade\tCluster_Seq\tCluster_Taxon\n')
        w.write('\n'.join(knn_summary_data))

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Preparing By-Catch Check Output Tables')

    df = pd.read_table(out_tsv)
    df = df[df.Major_Clade != 'No-clear-major-clade']
    euk_only = df[df['Domain'] == 'Eukaryota']

    df1 = df.groupby(['Query_Seq','Domain','Major_Group']).size().reset_index().rename(columns={0:'Frequency'})
    dmn_sprgrp = df1.groupby(['Domain','Major_Group']).size().reset_index().rename(columns={0:'Frequency'})
    dmn_sprgrp = dmn_sprgrp.sort_values(['Domain','Frequency'], ascending = [True, False])

    df2 = euk_only.groupby(['Query_Seq','Major_Group','Major_Clade','Minor_Clade']).size().reset_index().rename(columns={0:'Frequency'})
    euk_spr_mjr_mnr = df2.groupby(['Major_Group','Major_Clade','Minor_Clade']).size().reset_index().rename(columns={0:'Frequency'})
    euk_spr_mjr_mnr.rename(columns={'Major_Group': 'Super_Group'}, inplace = True)


    df3 = euk_only.groupby(['Query_Seq','Major_Group','Major_Clade']).size().reset_index().rename(columns={0:'Frequency'})
    euk_spr_mjr = df3.groupby(['Major_Group','Major_Clade']).size().reset_index().rename(columns={0:'Frequency'})
    euk_spr_mjr.rename(columns={'Major_Group': 'Super_Group'}, inplace = True)
    euk_spr_mjr = euk_spr_mjr.sort_values(['Super_Group','Frequency'], ascending = [True, False])

    by_taxon = df.groupby(['Query_Seq','Cluster_Taxon']).size().reset_index().rename(columns={0:'Frequency'})

    by_spr_grp = df.groupby(['Query_Seq','Major_Group']).size().reset_index().rename(columns={0:'Frequency'})
    by_spr_grp.rename(columns={'Major_Group': 'Super_Group'}, inplace = True)

    by_mjr = df.groupby(['Query_Seq','Major_Clade']).size().reset_index().rename(columns={0:'Frequency'})

    # Taxon should also report, domain, supergroup, major-clade
    common_cluster_taxa = by_taxon.Cluster_Taxon.value_counts()[:5]
    common_cluster_sprgrp = by_spr_grp.Super_Group.value_counts()[:5]
    common_cluster_mjr = by_mjr.Major_Clade.value_counts()[:5]

    out_dmn_sprgrp = f'{out_base_tsv}Domain_Supergroup.tsv'
    out_euk_spr_mjr_mnr = f'{out_base_tsv}Eukaryote.Summary.tsv'
    out_euk_spr_mjr = f'{out_base_tsv}Eukaryote.Supergroup_MajorClade.tsv'

    top_taxa = f'{out_base_tsv}Top_Clustered_Taxa.tsv'
    top_sprgrp = f'{out_base_tsv}Top_Clustered_Supergroups.tsv'
    top_mjr = f'{out_base_tsv}Top_Clustered_MajorClades.tsv'

    dmn_sprgrp.to_csv(out_dmn_sprgrp, sep = '\t', index = False)
    euk_spr_mjr_mnr.to_csv(out_euk_spr_mjr_mnr, sep = '\t', index = False)
    euk_spr_mjr.to_csv(out_euk_spr_mjr, sep = '\t', index = False)

    common_cluster_taxa.to_csv(top_taxa, sep = '\t')
    common_cluster_sprgrp.to_csv(top_sprgrp, sep = '\t')
    common_cluster_mjr.to_csv(top_mjr, sep = '\t')

    shutil.copy2(out_euk_spr_mjr, out_dir)
    shutil.copy2(out_dmn_sprgrp, out_dir)


def phylo_contam_eval_trees(
        start_time,
        contam_trees_summary: dict,
        query_taxon: str,
        phylo_dir: str,
        out_dir: str,
        blen_mode: str,
        verbose: bool = False):

    outdir = f'{phylo_dir}/Summary_Tables/'

    Path(outdir).mkdir(exist_ok = True, parents = True)

    out_base_tsv = f'{outdir}{query_taxon}.PhyG_PhyloBased_Contam.'

    out_contam_data = f'{out_base_tsv}tsv'

    with open(out_contam_data,'w+') as w:
        w.write('Taxon\tTaxon_Sequence\tGene_Family\tDomain\tMajor_Group\t' \
            'Major_Clade\tMinor_Clade\tNum_Sister_Taxa\tSister_Taxa\tSister_Sequences\tBranch_Length\t' \
            f'{blen_mode.title()}_Branch_Length\tBranch_Length_Evaluation\n')

        for i in contam_trees_summary:
            if i[0] == query_taxon:
                x = '\t'.join(i)
                w.write(f'{x}\n')

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Preparing By-Catch Check Output Tables')

    df = pd.read_table(out_contam_data)
    df = df[df['Taxon'] == query_taxon]

    df_short = df[df['Branch_Length_Evaluation'] == 'short']

    euk_only = df[df['Domain'] == 'Eukaryota']
    euk_only_short = euk_only[euk_only['Branch_Length_Evaluation'] == 'short']

    dmn_sprgrp = df.groupby(['Domain','Major_Group']).size().reset_index().rename(columns={0:'Frequency'})
    dmn_sprgrp = dmn_sprgrp.sort_values(['Domain','Frequency'], ascending = [True, False])

    dmn_sprgrp_short = df_short.groupby(['Domain','Major_Group']).size().reset_index().rename(columns={0:'Frequency'})
    dmn_sprgrp_short = dmn_sprgrp_short.sort_values(['Domain','Frequency'], ascending = [True, False])

    euk_spr_mjr = euk_only.groupby(['Major_Group','Major_Clade']).size().reset_index().rename(columns={0:'Frequency'})
    euk_spr_mjr.rename(columns={'Major_Group': 'Super_Group'}, inplace = True)
    euk_spr_mjr = euk_spr_mjr.sort_values(['Super_Group','Frequency'], ascending = [True, False])

    euk_spr_mjr_short = euk_only_short.groupby(['Major_Group','Major_Clade']).size().reset_index().rename(columns={0:'Frequency'})
    euk_spr_mjr_short.rename(columns={'Major_Group': 'Super_Group'}, inplace = True)
    euk_spr_mjr_short = euk_spr_mjr_short.sort_values(['Super_Group','Frequency'], ascending = [True, False])

    euk_spr_mjr_mnr = euk_only.groupby(['Major_Group','Major_Clade','Minor_Clade']).size().reset_index().rename(columns={0:'Frequency'})

    # euk_only.rename(columns={'Major_Group': 'Super_Group'}, inplace = True)

    # by_spr_grp = euk_only.Major_Group.value_counts()[:5]
    #
    # by_mjr = euk_only.Major_Clade.value_counts()[:5]

    by_sngl_sis = df[df.Sister_Taxa.str.count(';') < 1]
    by_sngl_sis = by_sngl_sis.groupby(['Domain','Major_Group','Major_Clade','Minor_Clade','Sister_Taxa']).size().reset_index().rename(columns={0:'Frequency'})
    by_sngl_sis_cnts = by_sngl_sis.value_counts()[:5]

    all_sister_taxa = [i for j in df.Sister_Taxa.to_list() for i in j.split(";") if len(set(j.split(";"))) == 1]
    common_single_sis = sorted(([i, all_sister_taxa.count(i)] for i in list(set(all_sister_taxa))), key = lambda x: -x[1])[:5]

    # out_dmn_sprgrp = f'{out_base_tsv}Domain_Supergroup.tsv'
    out_dmn_sprgrp_short = f'{out_base_tsv}Domain_Supergroup.Short_Branch_Length.tsv'

    out_euk_spr_mjr_mnr = f'{out_base_tsv}Eukaryote.Summary.tsv'
    # out_euk_spr_mjr = f'{out_base_tsv}Eukaryote.Supergroup_MajorClade.tsv'
    out_euk_spr_mjr_short = f'{out_base_tsv}Eukaryote.Supergroup_MajorClade.Short_Branch_Length.tsv'

    top_sprgrp = f'{out_base_tsv}Top_Supergroups.tsv'
    top_mjr = f'{out_base_tsv}Top_MajorClades.tsv'

    # out_all_cmn_sis = f'{out_base_tsv}Top_Common_Sister_Taxa.tsv'
    out_all_sngl_sis = f'{out_base_tsv}Top_Single_Sister_Taxa.tsv'

    # with open(out_all_cmn_sis, 'w+') as w:
    #     w.write('Sister_Taxon\tCount\n')
    #     for i in common_single_sis:
    #         w.write(f'{i[0]}\t{i[1]}\n')

    # dmn_sprgrp.to_csv(out_dmn_sprgrp, sep = '\t', index = False)
    dmn_sprgrp_short.to_csv(out_dmn_sprgrp_short, sep = '\t', index = False)

    euk_spr_mjr_mnr.to_csv(out_euk_spr_mjr_mnr, sep = '\t', index = False)
    # euk_spr_mjr.to_csv(out_euk_spr_mjr, sep = '\t', index = False)
    euk_spr_mjr_short.to_csv(out_euk_spr_mjr_short, sep = '\t', index = False)

    by_sngl_sis_cnts.to_csv(out_all_sngl_sis, sep = '\t')
    # by_spr_grp.to_csv(top_sprgrp, sep = '\t')
    # by_mjr.to_csv(top_mjr, sep = '\t')

    # shutil.copy2(out_dmn_sprgrp, out_dir)
    shutil.copy2(out_dmn_sprgrp_short, out_dir)

    # shutil.copy2(out_euk_spr_mjr, out_dir)
    shutil.copy2(out_euk_spr_mjr_short, out_dir)


def eval_contam(
    start_time,
    fasta_file: str,
    query_taxon: str,
    contam_dir: str,
    diag_db: str,
    gen_code: str = '1',
    taxon_code: str = None,
    transcripts: bool = True,
    genbank: bool = False,
    prokaryotic: bool = False,
    phylo_base: bool = False,
    ref_msa: str = None,
    ref_trees: str = None,
    blen_mode: str = 'median',
    clust_thresh: int = 3,
    eval_thresh: float = 0.6,
    threads: int = 4,
    verbose: bool = True) -> None:

    contam_dir, query_fas = prep_query_data(
                                start_time,
                                fasta_file,
                                query_taxon,
                                contam_dir,
                                diag_db,
                                gen_code,
                                taxon_code,
                                transcripts,
                                genbank,
                                prokaryotic,
                                threads,
                                verbose)

    if verbose:
        print('\n#--------- Contamination Evaluation --------#')

    if phylo_base:
        phylo_contam(
            start_time,
            contam_dir,
            query_fas,
            ref_msa,
            ref_trees,
            query_taxon,
            blen_mode = blen_mode,
            threads = threads,
            verbose = verbose)


    else:

        knn_contam_eval(
            start_time,
            contam_dir,
            query_fas,
            diag_db,
            query_taxon,
            clust_thresh = clust_thresh,
            eval_thresh = eval_thresh,
            verbose = verbose)

if __name__ == '__main__':
    try:
        fasta_file = sys.argv[1]
        query_taxon = sys.argv[2]
        diag_db = sys.argv[3]

    except:
        print('Usage:\n\n    phyg_contam.py [FASTA-FILE] [TAXON-NAME] [ORF-DB]\n\n')
        sys.exit(1)

    if 'Picozoa' in query_taxon:
        transcripts = False
    else:
        transcripts = True

    ref_msa = '/home/dr_x/Desktop/PhyG/PhyG_Data/PhyG_Contam_Data/PhyG_Contam_MSAs/'
    ref_tree = '/home/dr_x/Desktop/PhyG/PhyG_Data/PhyG_Contam_Data/PhyG_Contam_Trees/'
    phylo_base = False

    eval_contam(
        fasta_file,
        query_taxon,
        diag_db,
        transcripts = transcripts,
        phylo_base = phylo_base,
        ref_msa = ref_msa,
        ref_trees = ref_tree,
        threads = 24)
