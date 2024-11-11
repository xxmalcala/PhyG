#!/usr/bin/env python3

"""
Prepares and processes multi-sequence alignments with varying degrees of stringency.

Processing can include:
-- Size Distribution --
- using mean/median size distributions of peptides per gene family to remove
- peptides that are significantly shorter (33%) or longer (200%) the median OR
- mean length.

-- GUIDANCE2 (SLOW; Medium Stringency)--
- uses the third party tool, GUIDANCE2, with MAFFT to identify putative outlier
- sequences to be removed.

-- K-mer + N-gram Density (FAST; Medium Stringency)--
- Uses NLP approaches to infer k-mer n-gram frequencies, which are used for
- unsupervised outlier detection

-- Similarity Analyses (FAST; Low-HIGH Stringency [thresholds])--
- Uses NLP approaches to Infer distributions of "sentence" similarities for a
- given gene family to infer sets of "outlier" sequences to remove from a given
- UNALIGNED FASTA file, above a normalized threshold
- - Thresholds for scores are set to 3: this seems like a "good point"
- - Threshold guid:
- - Threshold = 2: stringent
- - Threshold = 4: moderate
- - Threshold = 10: permissive
-
- For small alignments, with low peptide diversity, this approach, with a
- stringen threshold, works very well!
-
- This approach can be coupled with ML-based outlier detection, which better
- matches K-mer and GUIDANCE2 based approaches and minimizes false outliers

-- Recommend AGAINST rescaling data if using k-mer based approaches. However, rescaling
-- is helpful if using similarity analyses + outlier detection


Outputs include alignments with "ALL" gaps, "divvied" alignments (if using Divvier),
and trimmed alignments (if desired).
"""

import glob, random, shutil, subprocess, sys

import numpy as np

from collections import defaultdict
from Levenshtein import distance
from pathlib import Path
from scipy.spatial.distance import pdist, squareform
from scipy.stats import iqr

from Bio import SeqIO


def gather_gfs(
        msa_dir: str,
        out_dir: str,
        prot_dir: str,
        taxon_list: str,
        gf_list: str,
        delim: str = '|',
        og_delim: str = '|') -> list:
    """
    Prepare FASTA files for a set of gene families and taxa.

    Parameters
    ----------
    out_dir:     output directory to store prepared gene families data
    prot_dir:    directory with stored peptide FASTA files with assigned gene families
    taxon_list:  list of taxa to include in the analysis
    gf_list:     list of gene families to include in the analysis
    og_delim:    string delimiter for OG-name to ease sequence grabbing

    Returns
    ----------
    init_gf_dir:  path to folder with data for query gene families and taxa
    gf_dict:      dictionary of all the initial gene families
    """
    gf_dict = defaultdict(list)

    init_gf_dir = f'{out_dir}_PhyG_Tree/Initial_Gene_Families/'

    prepped_gfs = 0

    Path(init_gf_dir).mkdir(exist_ok = True, parents = True)

    if msa_dir:
        for f in glob.glob(f'{msa_dir}/*.fa*'):
            og = f.rpartition("/")[-1].partition(".")[0]
            # Skip over sequences with translated internal stop codons
            gf_dict[og] += [i for i in SeqIO.parse(f,'fasta') if '*' not in i.seq]

    else:
        query_taxa = [i.rstrip() for i in open(taxon_list).readlines() if i != '']
        query_gfs = [i.rstrip() for i in open(gf_list).readlines() if i != '']


        # Checks all the FASTA files o- - - CHECK THE TUNING OF THIS USING OG5_126611f a given folder
        for f in glob.glob(f'{prot_dir}/*.fa*'):
            file_taxon = f.rpartition("/")[-1].partition(".")[0]

            # Check that these data are from a query taxon and query gene family
            if file_taxon in query_taxa:
                for seq in SeqIO.parse(f,'fasta'):

                    og = seq.id.rpartition(og_delim)[-1]
                    if og in query_gfs:
                        gf_dict[og].append(seq)

    # Ensures skipping over "empty" gene families
    for k, v in gf_dict.items():
        if len(v) > 0:
            prepped_gfs += 1
            SeqIO.write(v, f'{init_gf_dir}{k}.Initial.fasta','fasta')

    if prepped_gfs == 0:
        print('\nWARNING: No data were prepared. Ensure that the list of query ' \
            f'taxa are found in:\n   {prot_dir}\n')
        print('Also, ensure that the gene family names are ' \
            f'in the gene family list and are\nseparated by the last "{og_delim}" ' \
            'in your query data.\n')
        print('Exiting PhyG-Tree')
        sys.exit(1)

    return init_gf_dir, gf_dict


def run_mafft(
        out_dir: str,
        fasta_file: str,
        out_fasta: str,
        threads: int = 4) -> None:
    """
    Runs MAFFT alignment of single-gene MSAs.

    Parameters
    ----------
    out_dir:     output directory to store aligned FASTA file
    fasta_file:  Path to FASTA file to align with MAFFT
    out_fasta:   Path/filename to output aligned sequences
    threads:     number of CPU threads to use
    """

    mafft_cmd = f'mafft --quiet --thread {threads} --auto {fasta_file} > {out_dir}{out_fasta}'

    mafft_rslt = subprocess.Popen(mafft_cmd,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE,
                            universal_newlines = True,
                            shell = True)
    stdout, sderr = mafft_rslt.communicate()


def size_filter_genes(
        gf_filt_wd: str,
        gf_dict: dict,
        size_distro_median: bool,
        size_distro_min: float,
        size_distro_max: float) -> None:
    """
    Filters sequences from FASTA files that are too small/large.

    Parameters
    ----------
    gf_filt_wd:          Directory with FASTA files to filter
    gf_dict:             Dictionary of sequences for each query gene family
    size_distro_median:  Use the median sequence length per gene family (True; False = Mean)
    size_distro_min:     Minimum proportion of the median/mean seq length to filter sequences
    size_distro_max:    Maximum proportion of the median/mean seq length to filter sequences
    """

    size_distro_gfs = {}
    size_filt_gfs = 0

    for k, v in gf_dict.items():
        if len(v) < 4:
            continue

        slens = [len(i) for i in v]

        if size_distro_median:
            min_len = size_distro_min * np.median(slens)
            max_len = size_distro_max * np.median(slens)

        else:
            min_len = size_distro_min * np.mean(slens)
            max_len = size_distro_max * np.mean(slens)

        size_distro_gfs[k] = [i for i in v if min_len <= len(i) <= max_len]

    for k, v in size_distro_gfs.items():
        if len(v) > 3:
            SeqIO.write(v, f'{gf_filt_wd}/{k}.Size_Filt.fasta','fasta')

    size_distro_gfs.clear()


def ngram_from_seq(
        seq: str,
        ngram_len: int,
        overlap: bool) -> list:

    if overlap:
        return ' '.join([seq[n:n+ngram_len] for n in range(len(seq) - ngram_len + 1)])

    else:
        return ' '.join([seq[n:n+ngram_len] for n in range(0, len(seq) - ngram_len + 1, ngram_len)])


def prep_peptide_ngrams(
        fasta_file: str,
        ngram_len: int,
        overlap: bool) -> list:

    peptide_ngrams = {}

    for i in SeqIO.parse(fasta_file,'fasta'):
        if overlap:
            peptide_ngrams[i.id] = ngram_from_seq(f'{i.seq}', ngram_len, True)
        else:
            peptide_ngrams[i.id] = ngram_from_seq(f'{i.seq}', ngram_len, False)

    return peptide_ngrams


def calc_distance_matrix(seqs: list):
    transformed_seqs = np.array(seqs).reshape(-1,1)

    dist_matrix = pdist(transformed_seqs, lambda x,y: distance(x[0],y[0]))

    return squareform(dist_matrix)


def bootstrap_metrics(
        dist_matrix,
        nseqs: int = 100,
        reps: int = 100,
        threshold: int = 3) -> list:

    dm_median = [np.median(i) for i in dist_matrix]

    pseudo_reps = [(np.mean(i), np.std(i)) for i in np.random.choice(dm_median, size=(nseqs, len(dm_median)))]

    pseudo_mean = np.mean([i[0] for i in pseudo_reps])

    pseudo_std = np.mean([i[1] for i in pseudo_reps])

    bootstrap_outlier_seqs = [n for n in range(len(dm_median)) if threshold <= ((dm_median[n] - pseudo_mean)/pseudo_std)]

    outlier_seq_scores = [((dm_median[n] - pseudo_mean)/pseudo_std) for n in range(len(dm_median))]

    return bootstrap_outlier_seqs, outlier_seq_scores

def iqr_metrics(
        dist_matrix,
        threshold: int = 2) -> list:

    dm_median = [np.median(i) for i in dist_matrix]
    q1, q3 = np.percentile(dm_median, 25), np.percentile(dm_median, 75)

    iqr_outlier_seqs = []

    for n in range(len(dm_median)):
        scr_q1 = (q1-dm_median[n])/(q3-q1)
        scr_q3 = (q3-dm_median[n])/(q3-q1)

        if dm_median[n] < q1 and threshold > scr_q1:
            iqr_outlier_seqs.append(n)

        elif dm_median[n] > q1 and threshold < scr_q3:
            iqr_outlier_seqs.append(n)

    return iqr_outlier_seqs


def sen_sim_remove_outlier_seqs(
        fasta_file: str,
        bootstrap: bool = True,
        iqr: bool = False,
        nseqs: int = 100,
        reps: int = 100,
        threshold: int = 2,
        threads: int = 4) -> list:

# add option to couple with LOF



    clean_seqs, outlier_seqs = [], []

    all_seqs = [i for i in SeqIO.parse(fasta_file,'fasta')]

    seq_list = [f'{i.seq}' for i in all_seqs]

    dist_matrix = calc_distance_matrix(seq_list)

    if not iqr:
        bootstrap_outlier_seqs, outlier_seq_scores = bootstrap_metrics(
                                                        dist_matrix,
                                                        nseqs,
                                                        reps,
                                                        threshold)

        if bootstrap:
            for n in range(len(all_seqs)):
                if n in bootstrap_outlier_seqs:
                    outlier_seqs.append(all_seqs[n])
                else:
                    clean_seqs.append(all_seqs[n])
        else:
            oss = np.array(outlier_seq_scores).reshape(-1,1)

            lof_preds = cluster_lof(
                            oss,
                            True,
                            threads = 4)

            for n in range(len(all_seqs)):
                if lof_preds[n] == -1:
                    outlier_seqs.append(all_seqs[n])
                else:
                    clean_seqs.append(all_seqs[n])

    else:
        iqr_outlier_seqs = iqr_metrics(dist_matrix, threshold)

        for n in range(len(all_seqs)):
            if n in iqr_outlier_seqs:
                outlier_seqs.append(all_seqs[n])
            else:
                clean_seqs.append(all_seqs[n])

    return outlier_seqs, clean_seqs, outlier_seq_scores


def kmer_remove_outlier_seqs(
        fasta_file: str,
        ngram_len: int = 5,
        overlap: bool = True,
        threads: int = 4) -> list:

    from sklearn.feature_extraction.text import TfidfVectorizer

    clean_seqs, outlier_seqs = [], []

    all_seqs = [i for i in SeqIO.parse(fasta_file,'fasta')]

    peptide_ngrams = prep_peptide_ngrams(
                        fasta_file,
                        ngram_len,
                        overlap)

    tfidf = TfidfVectorizer(
                ngram_range = (1, 4),
                min_df = 5)

    peptide_tfidf = tfidf.fit_transform(peptide_ngrams.values()).toarray()

    lof_preds = cluster_lof(
                    peptide_tfidf,
                    False,
                    threads)

    for n in range(len(lof_preds)):
        if lof_preds[n] == -1:
            outlier_seqs.append(all_seqs[n])
        else:
            clean_seqs.append(all_seqs[n])

    return outlier_seqs, clean_seqs


def cluster_lof(
        query_data,
        rescale: bool = False,
        threads: int = 4):

    from sklearn.neighbors import LocalOutlierFactor

    neighbors = 20

    if rescale:
        from sklearn.preprocessing import RobustScaler

        query_data = RobustScaler().fit_transform(query_data)

    if len(query_data) < 20:
        neighbors = max(int(len(query_data)/2), 2)

    clf = LocalOutlierFactor(
            n_neighbors = neighbors,
            n_jobs = threads)

    preds = clf.fit_predict(query_data)

    # add option to save the clf

    return preds


def filter_families(
        out_dir: str,
        init_gf_dir: str,
        gf_dict: dict,
        msa_prog: str = 'MAFFT',
        size_distro: bool = True,
        size_distro_median: bool = True,
        size_distro_min: float = 0.33,
        size_distro_max: float = 2.00,
        guidance_filt: bool = False,
        guidance_filt_reps: int = 10,
        kmer_filt: bool = False,
        kmer_filt_len: int = 5,
        outlier_dist_filt: bool = False,
        outlier_dist_filt_boot: bool = True,
        outlier_dist_filt_iqr: bool = False,
        outlier_dist_filt_thresh: int = 2,
        outlier_dist_filt_nseqs: int = 100,
        outlier_dist_filt_reps: int = 100,
        threads = 4) -> str:
    """
    Prepare FASTA files for a set of gene families and taxa.

    Parameters
    ----------
    out_dir:     output directory to store prepared gene families data
    prot_dir:    directory with stored peptide FASTA files with assigned gene families
    taxon_list:  list of taxa to include in the analysis
    gf_list:     list of gene families to include in the analysis
    og_delim:    string delimiter for OG-name to ease sequence grabbing

    Returns
    ----------
    init_gf_dir:  path to folder with data for query gene families and taxa
    gf_dict:      dictionary of all the initial gene families
    """

    gf_filt_wd = init_gf_dir

    unaln_dir = f'{out_dir}_PhyG_Tree/Single_Gene_MSAs/UnAligned_MSAs/'

    if size_distro:
        gf_filt_wd = f'{out_dir}_PhyG_Tree/Filtering_Gene_Families/Size_Distribution/'

        Path(gf_filt_wd).mkdir(exist_ok = True, parents = True)

        size_filter_genes(
                gf_filt_wd,
                gf_dict,
                size_distro_median,
                size_distro_min,
                size_distro_max)

        gf_dict.clear()

    if guidance_filt:
        pass

    else:

        Path(unaln_dir).mkdir(exist_ok = True, parents = True)

        if (kmer_filt or guidance_filt or outlier_dist_filt):
            rm_seq_dir = f'{out_dir}_PhyG_Tree/Filtering_Gene_Families/Outlier_Seqs/'

            Path(rm_seq_dir).mkdir(exist_ok = True, parents = True)

            for fasta_file in glob.glob(f'{gf_filt_wd}/*.fa*'):
                out_fas = f'{fasta_file.rpartition("/")[-1].partition(".")[0]}'

                if outlier_dist_filt:
                    outlier_seqs, clean_seqs, seq_scores = sen_sim_remove_outlier_seqs(
                                                            fasta_file,
                                                            outlier_dist_filt_boot,
                                                            outlier_dist_filt_iqr,
                                                            outlier_dist_filt_nseqs,
                                                            outlier_dist_filt_reps,
                                                            outlier_dist_filt_thresh,
                                                            threads)

                elif kmer_filt:
                    print('KMERS!')
                    outlier_seqs, clean_seqs = kmer_remove_outlier_seqs(
                                                fasta_file,
                                                kmer_filt_len,
                                                True,
                                                threads)

                if clean_seqs:
                    clean_seqs_fasta = f'{out_fas}.Post_Filter.Clean.fasta'
                    SeqIO.write(clean_seqs, f'{unaln_dir}{clean_seqs_fasta}','fasta')

                if outlier_seqs:
                    removed_seqs_fasta = f'{out_fas}.Post_Filter.Removed.fasta'
                    SeqIO.write(outlier_seqs, f'{rm_seq_dir}{removed_seqs_fasta}','fasta')
        else:
            for fasta_file in glob.glob(f'{gf_filt_wd}/*.fa*'):
                out_fas = f'{fasta_file.rpartition("/")[-1].rpartition(".")[0]}.UnAligned.fasta'
                shutil.copy2(fasta_file, f'{unaln_dir}{out_fas}')


    return unaln_dir


        # remove all the "bad" sequences and make new FASTA files in a "filtered" folder
        # update `gf_filt_wd` as new filtered folder


    # May want a function that checks how gappy sequences are ... e.g., if they don't align over 50% or more of the alignment...
    # --> can be IQR too? Similar to OD-seq IQR scoring? but on the low end only?
    # -- --> might just want to keep an option to remove sequences that align across < 40 % of the alignment AFTER gap trimming
    # -- -- --> Remove from pre-alignment, re-align, re-gap trim



def align_msa_files(
        unaln_dir: str,
        msa_prog: str = 'MAFFT',
        threads: int = 4
        ):

    aln_dir = unaln_dir.replace('/UnAligned_MSAs/','/Aligned_MSAs/')
    Path(aln_dir).mkdir(exist_ok = True, parents = True)

    pre_align_fastas = glob.glob(f'{unaln_dir}*.fa*')

    if msa_prog == 'MAFFT':
        for paf in pre_align_fastas:
            out_fasta = f'{paf.rpartition("/")[-1].partition(".")[0]}.MAFFT_Aligned.fasta'
            run_mafft(
                    aln_dir,
                    paf,
                    out_fasta,
                    threads
                    )

        return aln_dir


def prep_gfs_eval_msas(
        out_dir: str,
        msa_dir: str,
        prot_dir: str,
        taxon_list: str,
        gf_list: str,
        delim: str = '|',
        og_delim: str = '|',
        msa_prog: str = 'MAFFT',
        size_distro: bool = True,
        size_distro_median: bool = True,
        size_distro_min: float = 0.33,
        size_distro_max: float = 2.00,
        guidance_filt: bool = False,
        guidance_filt_reps: int = 10,
        kmer_filt: bool = False,
        kmer_filt_len: int = 5,
        outlier_dist_filt: bool = True,
        outlier_dist_filt_boot: bool = True,
        outlier_dist_filt_iqr: bool = False,
        outlier_dist_filt_thresh: int = 3,
        outlier_dist_filt_nseqs: int = 100,
        outlier_dist_filt_reps: int = 100,
        threads: int = 4) -> str:

    init_gf_dir, gf_dict = gather_gfs(
                                msa_dir,
                                out_dir,
                                prot_dir,
                                taxon_list,
                                gf_list,
                                delim,
                                og_delim
                                )


    unaln_wd = filter_families(
                out_dir,
                init_gf_dir,
                gf_dict,
                msa_prog,
                size_distro,
                size_distro_median,
                size_distro_min,
                size_distro_max,
                guidance_filt,
                guidance_filt_reps,
                kmer_filt,
                kmer_filt_len,
                outlier_dist_filt,
                outlier_dist_filt_boot,
                outlier_dist_filt_iqr,
                outlier_dist_filt_thresh,
                outlier_dist_filt_nseqs,
                outlier_dist_filt_reps,
                threads)


    aln_wd = align_msa_files(
                unaln_wd,
                msa_prog,
                threads
                )


    return aln_wd


if __name__ == '__main__':
    try:
        gf_list = sys.argv[1]
        taxon_list = sys.argv[2]
    except:
        print('gf_list, taxa_list')
        sys.exit()

    out_dir = 'Should_Work'
    msa_dir = 'Should_Work_PhyG_Tree/TEST_Gene_Families/'
    # msa_dir = None
    prot_dir = '/home/dr_x/Desktop/PhyG/Prots_PhyG_Test'

    prep_gfs_eval_msas(
        out_dir,
        msa_dir,
        prot_dir,
        taxon_list,
        gf_list,
        og_delim = '|',
        msa_prog = 'MAFFT',
        size_distro = True,
        size_distro_median = True,
        size_distro_min = 0.33,
        size_distro_max = 2.00,
        guidance_filt = False,
        guidance_filt_reps = 10,
        kmer_filt = True,
        kmer_filt_len = 5,
        outlier_dist_filt = False,
        outlier_dist_filt_boot = True,
        outlier_dist_filt_iqr = False,
        outlier_dist_filt_thresh = 3,
        outlier_dist_filt_nseqs = 100,
        outlier_dist_filt_reps = 100,
        threads = 4)
