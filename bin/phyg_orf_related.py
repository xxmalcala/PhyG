#!/usr/bin/env python3

import subprocess

from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


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
    prots:          FASTA file contains protein sequences

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
        dmnd_cmd = dmnd_cmd.replace("blastx","blastp").rpartition(' qframe')[0]

    if min_id > 0:
        dmnd_cmd = dmnd_cmd.replace("sensitive ",f'sensitive --id {min_id} ')

    dmnd_rslt = subprocess.run(dmnd_cmd.split(),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True)

    return f'{out_dir}{out_tsv}'


def extract_orf(
        orf_coords: str,
        fasta_file: str,
        og_delim: str = '|',
        blastx: bool = True,
        blastp: bool = False,
        hmmer: bool = False) -> dict:
    """
    Extracts the ORF coordinates and returns FASTA files of the CDSs and corresponding
    translated CDSs.

    Parameters
    ----------
    orf_coords:  file with ORF coordinates
    fasta_file:  FASTA formatted file
    og_delim:    string delimiter for the hit -- the orthologous gene family ID needs to be after the final "split"!
    gen_code:    genetic code for translating CDSs
    blastx:      output from DIAMOND implementation of "BLASTX"
    hmmer:       output from hmmsearch

    Returns
    ----------
    orf_dict:  dictionary with extracted ORFs, their translations, and orthologous gene family
    """

    # Note to self (re: delimiter): Need an "Nth" delimiter option and option for last n-characters
    # Nth delimiter should be joined with "delim" as mandatory (make -1 by default)
    # n-characters would supercede delim/nth-delimiter

    if blastx:
        orf_dict = {line.split('\t')[0]:[(og_delim + line.split('\t')[1].rpartition(og_delim)[-1], int(line.split('\t')[-1]),
                    [int(line.split('\t')[4]), int(line.split('\t')[5])], float(line.split('\t')[-2]))] \
                    for line in open(orf_coords).readlines()}

        for i in SeqIO.parse(fasta_file, 'fasta'):
            if i.id in orf_dict:
                og, rf, coords, hsp = orf_dict[i.id][0]
                tmp_seq = i.seq[min(coords)-1:max(coords)]
                overhang = len(tmp_seq)%3
                if rf < 0:
                    tmp_seq = tmp_seq.reverse_complement()
                if overhang > 0:
                    tmp_seq = tmp_seq[:-overhang]
                orf_dict[i.id].append(tmp_seq)

    elif blastp:
        orf_dict = {line.split('\t')[0]:[(og_delim + line.split('\t')[1].rpartition(og_delim)[-1], float(line.split('\t')[-1]))] \
                    for line in open(orf_coords).readlines()}

        for i in SeqIO.parse(fasta_file, 'fasta'):
            if i.id in orf_dict:
                orf_dict[i.id].append(i.seq)

    return orf_dict


def translate_orfs(
        fasta_file: str,
        gen_code: str) -> str:

    translated_fasta = fasta_file.replace(".NTD.",".AA.")

    prot_seqs = []
    for i in SeqIO.parse(fasta_file, 'fasta'):
        i.seq = i.seq.translate(table = int(gen_code)).rstrip("*")
        prot_seqs.append(i)

    SeqIO.write(prot_seqs, translated_fasta, 'fasta')

    return translated_fasta


def finalize_orfs(
        out_dir: str,
        orf_dict: dict,
        taxon_name: str,
        delim: str = '|',
        top_hsp_only: bool = False,
        prots: bool = False) -> None:
    """
    Store outputs from OG-assignment and ORF-calling steps

    Parameters
    ----------
    out_dir:    directory to store FASTA files
    orf_dict:   dictionary with extracted ORFs, their translations, and orthologous gene family
    out_fasta:  FASTA filename to store the finalized ORFs/AAs
    delim:      delimiter for orthologous gene family assignment in sequence name
    prots:      orf_dict contains info for only protein sequences

    Returns
    ----------
    None
    """

    all_seqs = []

    top_hsp_seqs = []

    top_hits = []

    if top_hsp_only:
        d = defaultdict(list)
        for k, v in orf_dict.items():
            d[v[0][0]].append((k, v[0][-1]))

        for k, v in d.items():
            top_scr = max(v, key=lambda x: float(x[-1]))
            top_hits.append(top_scr[0])

    for k, v in orf_dict.items():
        seq_name = f'{k}{delim}{v[0][0]}'
        if top_hsp_only:
            if k in top_hits:
                top_hsp_seqs.append(SeqRecord(v[-1], id = seq_name, description = '', name = ''))

        all_seqs.append(SeqRecord(v[-1], id = seq_name, description = '', name = ''))

    if not prots:
        SeqIO.write(all_seqs, f'{out_dir}{taxon_name}.OG_Assigned.NTD.fasta', 'fasta')

        if top_hsp_only:
                SeqIO.write(top_hsp_seqs, f'{out_dir}{taxon_name}.OG_Assigned.Top_HSP.NTD.fasta', 'fasta')

                return f'{out_dir}{taxon_name}.OG_Assigned.Top_HSP.NTD.fasta'

        return f'{out_dir}{taxon_name}.OG_Assigned.NTD.fasta'

    else:
        SeqIO.write(all_seqs, f'{out_dir}{taxon_name}.OG_Assigned.AA.fasta', 'fasta')

        if top_hsp_only:
                SeqIO.write(top_hsp_seqs, f'{out_dir}{taxon_name}.OG_Assigned.Top_HSP.AA.fasta', 'fasta')

                return f'{out_dir}{taxon_name}.OG_Assigned.Top_HSP.AA.fasta'

        return f'{out_dir}{taxon_name}.OG_Assigned.AA.fasta'


def orf_call_diamond(
            out_dir: str,
            fasta_file: str,
            taxon_name: str,
            diamond_db: str,
            delim: str = '|',
            og_delim: str = '|',
            gen_code: str = '1',
            threads: int = 4,
            evalue: float = 1e-10,
            min_id: int = 0,
            max_hits: int = 1,
            min_hit_cover: int = 60,
            min_qry_cover: int = 0,
            top_hsp_only: bool = False,
            prots = False) -> str:

    out_tsv = f'{taxon_name}.DIAMOND_Hits.tsv'

    dmnd_tsv = diamond_og_assign(
                    out_dir,
                    out_tsv,
                    fasta_file,
                    diamond_db,
                    threads,
                    evalue,
                    min_id,
                    max_hits,
                    min_hit_cover,
                    min_qry_cover,
                    prots
                    )

    if not prots:
        orf_dict = extract_orf(
                        dmnd_tsv,
                        fasta_file,
                        og_delim
                        )

        if delim == og_delim:
            delim = ''

        out_ntds = finalize_orfs(
                    out_dir,
                    orf_dict,
                    taxon_name,
                    delim,
                    top_hsp_only,
                    prots
                    )

        return out_ntds, dmnd_tsv

    else:
        orf_dict = extract_orf(
                        dmnd_tsv,
                        fasta_file,
                        og_delim,
                        False,
                        True
                        )

        if delim == og_delim:
            delim = ''

        out_peps = finalize_orfs(
                    out_dir,
                    orf_dict,
                    taxon_name,
                    delim,
                    top_hsp_only,
                    True
                    )

        return out_peps, dmnd_tsv
