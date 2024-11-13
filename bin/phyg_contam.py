#!/usr/bin/env python3

import glob

from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from sklearn.feature_extraction.text import TfidfVectorizer

from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph

"""
Notes for self:

pull prep_peptide_ngrams and ngram_from_seq from phyg_msa_prep?

need taxonomy dictionary for the contamination database or update the taxon names
in the database and use NCBITaxa from ete3 to do the heavy lifting (might be more work)

Just does the basics for now... need to incorporate the 'one-offs' for the TIdeS
decontamination tests.

Still needs:
+ folder prep
+ ORF calling/handling
+ incorporate the contam-pipe functions into this script too

"""


def parse_og_fams(
        updated_og_dir: str):

    og_dict = defaultdict(list)


def prep_peptide_ngrams(
        gf_seqs: list,
        ngram_len: int = 5,
        overlap: bool = True) -> list:

    peptide_ngrams = {}

    for i in gf_seqs:
        if overlap:
            peptide_ngrams[i.id] = ngram_from_seq(f'{i.seq}', ngram_len, True)

        else:
            peptide_ngrams[i.id] = ngram_from_seq(f'{i.seq}', ngram_len, False)

    return peptide_ngrams


def ngram_from_seq(
        seq: str,
        ngram_len: int,
        overlap: bool) -> list:

    if overlap:
        return ' '.join([seq[n:n+ngram_len] for n in range(len(seq) - ngram_len + 1)])

    else:
        return ' '.join([seq[n:n+ngram_len] for n in range(0, len(seq) - ngram_len + 1, ngram_len)])


def knn_eval_contam(
        query_orf_fasta: str,
        og_reference_fasta: str,
        query_taxon: str):

    og_dict = defaultdict(list)

    for i in SeqIO.parse(query_orf_fasta,'fasta'):
        og_dict[i.id[-10:]].append(i)

    for i in SeqIO.parse(og_reference_fasta,'fasta'):
        if i.id[[-10:] in og_dict.keys():
            og_dict[i.id[-10:]].append(i)

    out_data = []

    tfidf = TfidfVectorizer(
                ngram_range = (1,5),
                min_df = 5)

    for k, v in og_dict.items():
        query_cluster_names = []
        eval_clusters = defaultdict(list)
        query_clusters = defaultdict(dict)

        knn_graph = kneighbors_graph(peptide_tfidf, 5, include_self=False)

        agg_mdl = AgglomerativeClustering(distance_threshold=2, n_clusters=None)

        peptide_ngrams = prep_peptide_ngrams(v)

        peptide_labels = list(peptide_ngrams.keys())

        peptide_tfidf = tfidf.fit_transform(peptide_ngrams.values()).toarray()

        agg_mdl = agg_mdl.fit(knn_graph.toarray())

        for n in range(len(model2.labels_)):
            s = peptide_labels[n]
            c = agg_mdl.labels_[n]

            eval_clusters[c].append(s)

            if s.startswith(query_taxon):
                query_cluster_names.append(c)

        for n in list(set(query_cluster_names)):
            s = eval_clusters[n]

            # For Katz-lab naming scheme! Needs to be updated!
            # Can use a simple dictionary if using a given database (save some time!)
            q_seqs = [i for i in s if query_taxon in i]
            c_seqs = list(set([i[:10] for i in s if query_taxon not in i]))

            if c_seqs:
                for q in q_seqs:
                    for i in c_seqs:
                        out_data.append(f'{q}\t{i[:2]}\t{i[:5]}\t{i}\n')

    with open('KNN_Test.Agglom_Thresh_2.TFIDF_1_5.Summary.tsv','w+') as w:
        w.write('Query_Sequence\tMajor_Clade\tMinor_Clade\tCluster_Taxon\n')
        w.write(''.join(out_data))
