#!/usr/bin/env python3

import subprocess

from pathlib import Path

from Bio import SeqIO


def remove_short_seqs(
        out_dir: str,
        taxon_name: str,
        fasta_file: str,
        min_length: int = 200,
        delim: str = '|',
        transcript: bool = True,
        wgs: bool = False) -> str:
    """
    Remove short sequences from FASTA file

    Parameters
    ----------
    out_dir:     output directory to store filtered data
    taxon_name:  species/taxon name or abbreviated code
    fasta_file:  FASTA formatted file
    min_len:     minimum ORF length to consider [default is 200nt]
    delimiter:   string to separate portions of the sequence name [default is "|"]
    transcript:

    Returns
    ----------
    size_filt_fas:  FASTA formatted file with short sequences removed
    """
    seq_num = 1

    size_filt_seqs = []
    name_conversion = {}

    size_filt_fas = f'{out_dir}{taxon_name}.{min_length}bp.fasta'

    for i in SeqIO.parse(fasta_file,'fasta'):
        if len(i) >= min_length:
            seq_name = f'{taxon_name}{delim}'
            if wgs:
                # Append accession information for the ORF -- to do
                pass
            elif transcript:
                seq_name += f'Transcript_{seq_num}_Length_{len(i)}'
                seq_num += 1
                if '_cov_' in i.id:
                    kcov = f'{float(i.id.partition("_cov_")[-1].split("_")[0]):.2f}'
                    seq_name += f'_Cov_{kcov}'
            name_conversion[i.id] = seq_name
            i.id = seq_name
            i.description = ''
            i.name = ''
            size_filt_seqs.append(i)

    SeqIO.write(size_filt_seqs, size_filt_fas, 'fasta')

    return size_filt_fas


def run_barrnap(fasta_file: str, threads: int) -> list:
    """
    Run Barrnap to remove easily identifiable rRNAs.

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    threads:     number of cpu threads to use

    Returns
    ----------
    rRNA_seqs:  list of putative rRNA sequence names
    """

    rRNA_seqs = []
    kdm = ['bac','arc','mito','euk']

    for k in kdm:
        bnp_cmd = ['barrnap', '--kingdom', k, '--threads', f'{threads}', fasta_file]

        bnp_rslt = subprocess.run(
                        bnp_cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True
                        )

        rRNA_seqs += [i.split('\t')[0] for i in bnp_rslt.stdout.split('\n') if '|' in i]

    return rRNA_seqs


def remove_rRNA(out_dir: str, taxon_name: str, fasta_file: str,  min_len: int, threads: int) -> str:
    """
    Remove putative rRNA sequences from FASTA file

    Parameters
    ----------
    out_dir:     output directory to store filtered data
    taxon_name:  species/taxon name or abbreviated code
    fasta_file:  FASTA formatted file
    min_len:     minimum ORF length to consider [default is 300nt]
    threads:     number of cpu threads to use

    Returns
    ----------
    rRNA_clean_fas:  FASTA formatted file without putative rRNA sequences
    """

    rRNA_seqs = run_barrnap(fasta_file, threads)

    rRNA_clean_fas = f'{out_dir}{taxon_code}.{min_len}bp.Filt_rRNA.fas'
    rRNA_fas = f'{out_dir}{taxon_code}.{min_len}bp.rRNA_Seqs.fas'

    rRNA_contam, clean_seqs = [],[]

    # bin sequences as either putative rRNAs or "clean"
    for i in SeqIO.parse(fasta_file, 'fasta'):
        if i.id in rRNA_seqs:
            rRNA_contam.append(i)
        else:
            clean_seqs.append(i)

    SeqIO.write(rRNA_contam, rRNA_fas, 'fasta')
    SeqIO.write(clean_seqs, rRNA_clean_fas, 'fasta')

    return rRNA_clean_fas
