#!/usr/bin/env python3

import subprocess, time
import numpy as np

import sys

from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from bin import phyg_gen_util as gut

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
        prots: bool = False,
        eval_gen_code: bool = False) -> str:
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

    if eval_gen_code:
        dmnd_cmd  = f'diamond blastx --very-sensitive ' \
                    f'-q {fasta_file} ' \
                    f'-d {diamond_db} ' \
                    f'-p {threads} ' \
                    f'-e {evalue} ' \
                    f'-k {max_hits} ' \
                    f'--subject-cover {min_hit_cover} ' \
                    f'-o {out_tsv} ' \
                    f'-f 5'

    if gut.check_file_exists(f'{out_dir}{out_tsv}'):
        return f'{out_dir}{out_tsv}'

    dmnd_rslt = subprocess.run(dmnd_cmd.split(),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True)
    if not eval_gen_code:
        return f'{out_dir}{out_tsv}'


def hmmer_og_assign(
        out_dir: str,
        out_tsv: str,
        fasta_file: str,
        hmmer_db: str,
        threads: int = 4,
        evalue: float = 1e-10) -> str:
    """
    Assign gene families with Hmmsearch against an hmm-database.

    Parameters
    ----------
    out_dir:     directory to store output TSV file of OG hits
    out_tsv:     output filenames
    fasta_file:  FASTA formatted file
    hmmer_db:    path to DIAMOND formatted database
    threads:     number of threads allowed for DIAMOND

    Returns
    ----------
    out_tsv:  path to TSV file with reported hits
    """

    # Worth noting that hmmsearch is I/O bound, so more threads than 12
    # does not lead to significant speed-ups

    hmmer_cmd  = f'hmmsearch ' \
                f'-E {evalue} ' \
                f'--cpu {min(threads, 12)} ' \
                f'--tblout {out_dir}{out_tsv} ' \
                f'-A {out_dir}{out_tsv.replace(".tsv",".sto")} ' \
                f'{hmmer_db} ' \
                f'{fasta_file}'

    if not gut.check_file_exists(f'{out_dir}{out_tsv.replace(".tsv",".sto")}'):

        hmmer_rslt = subprocess.run(hmmer_cmd.split(),
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True)

    return f'{out_dir}{out_tsv}', f'{out_dir}{out_tsv.replace(".tsv",".sto")}'


def parse_hmmer_search(
        hmmer_tblout: str,
        hmmer_sto: str,
        delim: str = '|',
        prots: bool = False) -> dict:
    """
    Parses the hmmer outputs, returning a dictionary of relevant information.

    Parameters
    ----------
    hmmer_tblout:  tblout from hmmsearch
    hmmer_sto:     stockholm format alignments from hmmsearch

    Returns
    ----------
    hmmer_hits:  summary of the hmmsearch hits, includes alignment coordinates
    """

    hmmer_hits = defaultdict(list)

    hmmer_coords = {i.split()[1].split("/")[0]:i.split()[1].split("/")[1] for i in open(hmmer_sto).readlines() if i.startswith('#=GS ')}

    for line in open(hmmer_tblout).readlines()[3:-10]:
        seq_info = list(np.array(line.split())[[0,2,4,5]])

        seq_info += [int(i)*3 for i in hmmer_coords[seq_info[0]].split('-')]

        if not prots:
            hmmer_hits[line.split()[0].rpartition(delim)[0]].append(seq_info)
        else:
            hmmer_hits[line.split()[0]].append(seq_info)

    return hmmer_hits


def extract_orf(
        orf_coords: str,
        fasta_file: str,
        og_delim: str = '|',
        blastx: bool = True,
        blastp: bool = False,
        hmmer_tblout: str = None,
        hmmer: bool = False,
        prots: bool = False) -> dict:
    """
    Extracts the ORF coordinates and returns FASTA files of corresponding CDSs

    Parameters
    ----------
    orf_coords:  file with ORF coordinates
    fasta_file:  FASTA formatted file
    og_delim:    string delimiter for the hit -- the orthologous gene family ID needs to be after the final "split"!
    gen_code:    genetic code for translating CDSs
    blastx:        output from DIAMOND implementation of "BLASTX"
    hmmer_tblout:  tblout from hmmsearch
    hmmer_sto:     stockholm format alignments from hmmsearch
    hmmer:       outputs from hmmsearch

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

    elif hmmer:
        orf_dict = {}
        hmmer_dict = parse_hmmer_search(
                    hmmer_tblout,
                    orf_coords,
                    og_delim,
                    prots)

        for i in SeqIO.parse(fasta_file,'fasta'):
            if i.id in hmmer_dict:
                top_hit = max(hmmer_dict[i.id], key = lambda x: float(x[3]))
                qstart, qend = top_hit[-2:]

                if not prots:
                    rf = int(top_hit[0].rpartition('RF')[-1])
                    if rf < 4:
                        qstart -= (4 - rf)
                        qend += (rf - 1)
                        qorf = i.seq[qstart:qend]
                    else:
                        qstart -= (7 - rf)
                        qend += (rf - 4)
                        qorf = i.seq.reverse_complement()[qstart:qend]

                else:
                    qorf = i.seq
                    rf = 1

                if qorf:
                    i.seq = qorf
                    i.description = ''
                    i.name = ''
                    orf_dict[i.id] = [(f'{og_delim}{top_hit[1]}', rf, [qstart, qend], float(top_hit[3])), i.seq]

    return orf_dict


def translate_all_frames(
        seq,
        gen_code: str = '1') -> list:

    seq_all_frames = []
    for n in range(1,7):
        seq_name = f'{seq.id}|RF{n}'
        if n < 4:
            if len(seq.seq[n-1:])%3 > 0:
                rf_seq = seq.seq[n-1:-(len(seq.seq[n-1:])%3)].translate(int(gen_code))
            else:
                rf_seq = seq.seq[n-1:].translate(int(gen_code))

        else:
            if len(seq.seq[n-4:])%3 > 0:
                rf_seq = seq.seq.reverse_complement()[n-4:-(len(seq.seq[n-4:])%3)].translate(table = int(gen_code))
            else:
                rf_seq = seq.seq.reverse_complement()[n-4:].translate(table = int(gen_code))

        seq_all_frames.append(f'>{seq_name}\n{rf_seq}')

    return seq_all_frames


def translate_orfs(
        fasta_file: str,
        gen_code: str,
        all_frames_fasta = None,
        hmmer: bool = False) -> str:

    if not hmmer:
        translated_fasta = fasta_file.replace(".NTD.",".AA.")

        if gut.check_file_exists(translated_fasta):
            return translated_fasta

        prot_seqs = []
        for i in SeqIO.parse(fasta_file, 'fasta'):
            i.seq = i.seq.translate(table = int(gen_code)).rstrip("*")
            prot_seqs.append(i)

        SeqIO.write(prot_seqs, translated_fasta, 'fasta')

        return translated_fasta

    else:
        seqs_all_frames = []

        for i in SeqIO.parse(fasta_file,'fasta'):
            seqs_all_frames += translate_all_frames(i, gen_code)

        with open(all_frames_fasta,'w+') as w:
            w.write('\n'.join(seqs_all_frames))


def finalize_orfs(
        out_dir: str,
        orf_dict: dict,
        taxon_name: str,
        delim: str = '|',
        top_hsp_only: bool = False,
        prots: bool = False,
        hmmer: bool = False) -> None:
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
        start_time,
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
        prots: bool = False,
        verbose: bool = False) -> str:

    out_tsv = f'{taxon_name}.DIAMOND_Hits.tsv'

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Assigning Gene-Families')

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
                    prots)

    if not prots:
        if verbose:
            print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Extracting ORFs')

        orf_dict = extract_orf(
                        dmnd_tsv,
                        fasta_file,
                        og_delim)

        if delim == og_delim:
            delim = ''

        if verbose:
            print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Finalizing Initial ORFs')

        out_ntds = finalize_orfs(
                    out_dir,
                    orf_dict,
                    taxon_name,
                    delim,
                    top_hsp_only,
                    prots)

        return out_ntds, dmnd_tsv

    else:
        orf_dict = extract_orf(
                        dmnd_tsv,
                        fasta_file,
                        og_delim,
                        False,
                        True)

        if delim == og_delim:
            delim = ''

        if verbose:
            print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Finalizing Assignments')

        out_peps = finalize_orfs(
                    out_dir,
                    orf_dict,
                    taxon_name,
                    delim,
                    top_hsp_only,
                    True)

        return out_peps, dmnd_tsv


def orf_call_hmmer(
        start_time,
        out_dir: str,
        fasta_file: str,
        taxon_name: str,
        hmmer_db: str,
        delim: str = '|',
        og_delim: str = '|',
        gen_code: str = '1',
        threads: int = 4,
        evalue: float = 1e-10,
        top_hsp_only: bool = False,
        prots: bool = False,
        verbose: bool = False) -> str:

    out_tsv = f'{taxon_name}.HMMER_Hits.tsv'

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Tanslating in All Reading Frames')

    if not prots:
        peptide_fasta = f'{out_dir}{taxon_name}.All_Six_Frames.fasta'

        if not gut.check_file_exists(peptide_fasta):
            translate_orfs(
                fasta_file,
                gen_code,
                peptide_fasta,
                hmmer = True)

    else:
        peptide_fasta = fasta_file

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Assigning Gene-Families')

    hmmer_tblout, hmmer_sto = hmmer_og_assign(
                                out_dir,
                                out_tsv,
                                peptide_fasta,
                                hmmer_db,
                                threads,
                                evalue)
    # if not prots:
    #     Path.unlink(Path.cwd() / peptide_fasta)

    if not prots:
        if verbose:
            print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Extracting ORFs')

        orf_dict = extract_orf(
                        hmmer_sto,
                        fasta_file,
                        og_delim,
                        hmmer_tblout = hmmer_tblout,
                        hmmer = True,
                        blastx = False)

        if delim == og_delim:
            delim = ''

        if verbose:
            print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Finalizing Initial ORFs')

        out_ntds = finalize_orfs(
                    out_dir,
                    orf_dict,
                    taxon_name,
                    delim,
                    top_hsp_only,
                    prots)

        return out_ntds

    else:
        orf_dict = extract_orf(
                        hmmer_sto,
                        fasta_file,
                        og_delim,
                        hmmer_tblout = hmmer_tblout,
                        hmmer = True,
                        blastx = False,
                        prots = True)

        if delim == og_delim:
            delim = ''

        if verbose:
            print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Finalizing Assignments')

        out_peps = finalize_orfs(
                    out_dir,
                    orf_dict,
                    taxon_name,
                    delim,
                    top_hsp_only,
                    True)

        return out_peps
