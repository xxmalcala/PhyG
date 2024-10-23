#!/usr/bin/env python3

"""
Script is intended to capture the longest isoform from a given FASTA file
harboring nucleotide CDSs downloaded from RefSeq/GenBank.
"""


import sys

from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def assess_seq_info(seq_description: str) -> list:
    """
    Returns locus information for a given sequence, based on its "description"

    Parameters
    ----------
    seq_description:  full descriptive name, presumably from GenBank

    Returns
    ----------
    loc_info:   locus for the given protein-coding gene
    prot_info:  unique accession code for the given protein-coding gene
    """

    if 'protein_id=' not in seq_description:
        return None, None

    elif 'locus_tag=' in seq_description:
        loc_info = seq_description.partition('locus_tag=')[2].partition(']')[0]
        prot_info = seq_description.partition('protein_id=')[2].partition(']')[0]
        return loc_info, prot_info

    elif 'gene=' in seq_description:
        loc_info = seq_description.partition('gene=')[2].partition(']')[0]
        prot_info = seq_description.partition('protein_id=')[2].partition(']')[0]
        return loc_info, prot_info

    else:
        prot_info = seq_description.partition('protein_id=')[2].partition(']')[0]
        return None, prot_info



def check_complete_seq(seq: str, prok: bool = False) -> bool:
    """
    Checks that the given sequence is complete with start and termination codons

    Parameters
    ----------
    seq:   CDS sequence to evaluate
    prok:  whether or not to check if the sequence starts with major prokaryotic alternative codons

    Returns
    ----------
    True if sequence is a complete ORF, with start and stop codons, otherwise False
    """

    if len(seq)%3 == 0 and seq[-3:].upper() in ['TGA','TAA','TAG']:
        if not prok and seq[:3].upper() == 'ATG':
            return True

        # Prokaryotes can use a couple other initiation codons (TTG/GTG are most common)
        elif prok and seq[:3].upper() in ['ATG', 'TTG','GTG']:
            return True

        else:
            return False

    else:
        return False


def capture_isoforms(fasta_file: str, prok: bool = False) -> dict:
    """
    Captures isoforms for each locus from a GenBank-sourced FASTA file of CDSs

    Parameters
    ----------
    fasta_file:   FASTA file of GenBank/RefSeq CDS sequence to evaluate
    prok:         whether or not to check if the sequence starts with major prokaryotic alternative codons

    Returns
    ----------
    iso_dict:    dictionary of isoforms for each unique locus
    """

    iso_dict = defaultdict(list)
    s = 0
    cs = 0

    for i in open(fasta_file).read().split('>lcl|')[1:]:
        s += 1
        qseq = i.rpartition(']')[2].replace('\n','')

        if check_complete_seq(qseq, prok):
            cs += 1
            loc_info, prot_info = assess_seq_info(i.partition('\n')[0])

            if loc_info:
                iso_dict[loc_info].append((prot_info, qseq))

            elif prot_info:
                iso_dict[loc_info].append((prot_info, qseq))

    print(f'Found {s} CDSs, {cs} are complete, from {len(iso_dict)} unique loci')

    return iso_dict


def save_isoforms(fasta_file: str, longest_isoform: bool = True, prok: bool = False) -> list:
    """
    Keep the most useful isoforms based on user preference (either all of them, or the
    longest isoform for a given locus)

    Parameters
    ----------
    fasta_file:       FASTA file of GenBank/RefSeq CDS sequence to evaluate
    longest_isoform:  whether to keep just the longest isoform (True), otherwise all isoforms (False)
    prok:             whether or not to check if the sequence starts with major prokaryotic alternative codons

    Returns
    ----------
    all_isoforms:    list of all the isoforms to keep for further analyses
    """

    all_isoforms = []

    iso_dict = capture_isoforms(fasta_file, prok)



    for k, v in iso_dict.items():
        if longest_isoform:
            prot_id, seq_str = max(v, key=lambda x: len(x[1]))

            all_isoforms.append(SeqRecord(Seq(seq_str), prot_id, '', ''))
        else:
            for i in v:
                all_isoforms.append(SeqRecord(Seq(i[1]), i[0], '', ''))

    return all_isoforms
