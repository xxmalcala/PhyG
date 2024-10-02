#!/usr/bin/env python3

import subprocess

from Bio import SeqIO


def diamond_og_assign(
        out_dir: str,
        out_tsv: str,
        fasta_file: str,
        diamond_db: str,
        threads: int = 4,
        evalue: float = 1e-10,
        min_id: int = 0,
        max_hits: int = 1,
        min_hit_cover: int = 60,
        min_qry_cover: int = 0,
        prots = False) -> str:
    """
    DIAMOND BLAST(X/P) against a user-provided OG database

    Parameters
    ----------
    out_dir:        directory to store output TSV file of OG hits
    out_tsv:        output TSV filename
    fasta_file:     FASTA formatted file
    diamond_db:     path to DIAMOND formatted database
    threads:        number of threads allowed for DIAMOND
    evalue:         maximum e-value to report a "hit"
    min_id:         minimum percent identity to report a "hit"
    max_hits:       maximum number of matches to report per query sequence
    min_hit_cover:  minimum proportion of the subject covered by the query sequence
    min_qry_cover:  minimum proportion of the query sequence covered

    Returns
    ----------
    out_tsv:  path to TSV file with reported hits
    """

    dmnd_cmd  = f'diamond blastx --very-sensitive ' \
                f'-q {fasta_file} ' \
                f'-d {diamond_db} ' \
                f'-p {threads} ' \
                f'-e {evalue} ' \
                f'-k {max_hits} ' \
                f'--subject-cover {min_hit_cover} ' \
                f'--query-cover {min_qry_cover} ' \
                f'-o {out_dir}{out_tsv} ' \
                f'-f 6 qseqid sseqid length pident qstart qend sstart send evalue bitscore qframe'

    if prots:
        dmnd_cmd = dmnd_cmd.replace("blastx","blastp").rstrip(' qframe')

    if min_id > 0:
        dmnd_cmd = dmnd_cmd.replace("sensitive ",f'sensitive --id {min_id} ')

    dmnd_rslt = subprocess.run(dmnd_cmd.split(),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True)

    return out_tsv


def extract_orf(
        orf_coords: str,
        fasta_file: str,
        delim: str = '|',
        gen_code: str = '1',
        blastx: bool = True,
        hmmer: bool = False) -> dict:
    """
    Extracts the ORF coordinates and returns FASTA files of the CDSs and corresponding
    translated CDSs.

    Parameters
    ----------
    orf_coords:  file with ORF coordinates
    fasta_file:  FASTA formatted file
    delim:       string delimiter for the hit -- the orthologous gene family ID needs to be after the final "split"!
    gen_code:    genetic code for translating CDSs
    blastx:      output from DIAMOND implementation of "BLASTX"
    hmmer:       output from hmmsearch

    Returns
    ----------
    orf_dict:  dictionary with extracted ORF and translation
    """

    if blastx:
        orf_dict = {line.split('\t')[0]:[(int(line.split('\t')[-1]), line.split('\t')[1].rpartition(delim)[-1], [int(line.split('\t')[4]), int(line.split('\t')[5])])] \
                for line in open(orf_coords).readlines()}

        for i in SeqIO.parse(fasta_file, 'fasta'):
            if i.id in orf_dict:
                rf, og, coords = orf_dict[i.id][0]
                tmp_seq = i.seq[min(coords)-1:max(coords)]
                overhang = len(tmp_seq)%3
                if rf < 0:
                    tmp_seq = tmp_seq.reverse_complement()
                if overhang > 0:
                    tmp_seq = tmp_seq[:-overhang]
                orf_dict[i.id].append(tmp_seq)
                orf_dict[i.id].append(tmp_seq.translate(gen_code))

    return orf_dict
