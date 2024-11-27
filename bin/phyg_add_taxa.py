#!/usr/bin/env python3

import logging, shutil, subprocess, sys, time

from datetime import datetime, timedelta
from collections import defaultdict

from pathlib import Path

# from bin import... MAKE SURE TO ADD THIS AFTER FINISHED TESTING!!!
from bin import phyg_orf_related as oc

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



def check_complete_seq(
        seq: str,
        prok: bool = False) -> bool:
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


def capture_isoforms(
        start_time,
        fasta_file: str,
        prok: bool = False,
        verbose: bool = False) -> dict:
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

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Found {s} CDSs, {cs} are complete, from {len(iso_dict)} unique loci')

    return iso_dict


def save_isoforms(
        start_time,
        fasta_file: str,
        longest_isoform: bool = True,
        prok: bool = False,
        verbose: bool = False) -> list:
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

    iso_dict = capture_isoforms(
                start_time,
                fasta_file,
                prok,
                verbose)

    for k, v in iso_dict.items():
        if longest_isoform:
            prot_id, seq_str = max(v, key=lambda x: len(x[1]))

            all_isoforms.append(SeqRecord(Seq(seq_str), prot_id, '', ''))
        else:
            for i in v:
                all_isoforms.append(SeqRecord(Seq(i[1]), i[0], '', ''))

    return all_isoforms


def back_up_file(
        out_dir: str,
        fasta_file: str) -> None:

    Path(out_dir).mkdir(exist_ok = True, parents = True)

    shutil.copy2(fasta_file, out_dir)


def eval_wgs_genbank(
        start_time,
        out_dir: str,
        taxon_name: str,
        fasta_file: str,
        taxon_code: str = None,
        prokaryotic: bool = False,
        remove_isoforms: bool = True,
        delim: str = '|',
        verbose: bool = False) -> str:

    Path(out_dir).mkdir(exist_ok = True, parents = True)

    final_filt_fas = f'{out_dir}{taxon_name}.LongestIsoform'

    if taxon_code:
        final_filt_fas = final_filt_fas.replace(taxon_name, taxon_code)

    if not remove_isoforms:
        final_filt_fas.replace(".LongestIsoform",".ORFs")

    final_seqs = []
    name_conversion = {}

    if remove_isoforms:
        isoforms_to_keep = save_isoforms(
                            start_time,
                            fasta_file,
                            True,
                            prokaryotic,
                            verbose)

        if not isoforms_to_keep:
            print('ERROR: Unable to identify NCBI loci and remove redundant isoforms'
            f' in {fasta_file.rpartition("/")[-1]}. Please ensure that this file was sourced'
            ' from GenBank/RefSeq.')

            print("If this warning persists, please open an issue on PhyG's GitHub page.")
            sys.exit(1)

    else:
        isoforms_to_keep = save_isoforms(
                                start_time,
                                fasta_file,
                                False,
                                prokaryotic,
                                verbose)

    for i in isoforms_to_keep:
        seq_name = f'{taxon_name}{delim}{i.id}'
        if taxon_code:
            seq_name = f'{taxon_code}{delim}{seq_name}'

        name_conversion[i.id] = seq_name
        i.id = seq_name
        i.description = ''
        i.name = ''

        final_seqs.append(i)

    SeqIO.write(final_seqs, f'{final_filt_fas}.NTD.fasta', 'fasta')

    with open(f'{final_filt_fas}.NameConversion.tsv', 'w+') as w:
        w.write('Initial_ID\tUpdated_Name\n')
        for k, v in name_conversion.items():
            w.write(f'{k}\t{v}\n')

    return f'{final_filt_fas}.NTD.fasta'


def eval_wgs_non_gbk(
        out_dir: str,
        taxon_name: str,
        fasta_file: str,
        taxon_code: str,
        prokaryotic: bool = False,
        delim: str = '|') -> str:

    Path(out_dir).mkdir(exist_ok = True, parents = True)

    final_filt_fas = f'{out_dir}{taxon_name}.ORFs'

    if taxon_code:
        final_filt_fas = final_filt_fas.replace(taxon_name, taxon_code)

    final_seqs = []
    name_conversion = {}

    for i in SeqIO.parse(fasta_file,'fasta'):
        if len(i.seq)%3 == 0:
            seq_name = f'{taxon_name}{delim}{i.id.split(delim)[-1]}'

            if taxon_code:
                seq_name = f'{taxon_code}{delim}{seq_name}'

            name_conversion[i.id] = seq_name
            i.id = seq_name
            i.description = ''
            i.name = ''
            final_seqs.append(i)

        else:
            pass

    SeqIO.write(final_seqs, f'{final_filt_fas}.NTD.fasta', 'fasta')

    with open(f'{final_filt_fas}.NameConversion.tsv', 'w+') as w:
        w.write('Initial_ID\tUpdated_Name\n')
        for k, v in name_conversion.items():
            w.write(f'{k}\t{v}\n')

    return f'{final_filt_fas}.NTD.fasta'


def remove_short_seqs(
        out_dir: str,
        taxon_name: str,
        taxon_code: str,
        fasta_file: str,
        min_length: int = 200,
        delim: str = '|') -> str:
    """
    Remove short sequences from FASTA file

    Parameters
    ----------
    out_dir:     output directory to store filtered data
    taxon_name:  species/taxon name or abbreviated code
    fasta_file:  FASTA formatted file
    min_len:     minimum ORF length to consider [default is 200nt]
    delimiter:   string to separate portions of the sequence name [default is "|"]

    Returns
    ----------
    size_filt_fas:  FASTA formatted file with short sequences removed
    """

    Path(out_dir).mkdir(exist_ok = True, parents = True)

    seq_num = 1

    size_filt_seqs = []
    name_conversion = {}

    size_filt_fas = f'{out_dir}{taxon_name}.{min_length}bp.fasta'
    if taxon_code:
        size_filt_fas = size_filt_fas.replace(f'{out_dir}{taxon_name}',f'{out_dir}{taxon_code}')

    for i in SeqIO.parse(fasta_file,'fasta'):
        if len(i) >= min_length:
            seq_name = f'{taxon_name}{delim}'

            if taxon_code:
                seq_name = f'{taxon_code}{delim}{seq_name}'

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

    with open(size_filt_fas.replace(".fasta",".NameConversion.tsv"), 'w+') as w:
        w.write('Initial_ID\tUpdated_Name\n')
        for k, v in name_conversion.items():
            w.write(f'{k}\t{v}\n')

    return size_filt_fas


def run_barrnap(
        fasta_file: str,
        threads: int) -> list:
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


def remove_rRNA(
        out_dir: str,
        taxon_name: str,
        fasta_file: str,
        min_len: int,
        threads: int) -> str:
    """
    Remove putative rRNA sequences from FASTA file

    Parameters
    ----------
    out_dir:     output directory to store filtered data
    taxon_name:  taxon/sample name or abbreviated code
    fasta_file:  FASTA formatted file
    min_len:     minimum ORF length to consider [default is 300nt]
    threads:     number of cpu threads to use

    Returns
    ----------
    rRNA_clean_fas:  FASTA formatted file without putative rRNA sequences
    """

    Path(out_dir).mkdir(exist_ok = True, parents = True)

    rRNA_seqs = run_barrnap(fasta_file, threads)

    rRNA_clean_fas = f'{out_dir}{taxon_name}.{min_len}bp.rRNA_Filtered.fasta'
    rRNA_fas = f'{out_dir}{taxon_name}.{min_len}bp.rRNA_Seqs.fasta'

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


def collapse_isoforms(
        out_dir: str,
        fasta_file: str,
        pid: float = 0.97,
        threads: int = 4) -> str:

    Path(out_dir).mkdir(exist_ok = True, parents = True)

    cdhit_out_fasta = f'{out_dir}{fasta_file.rpartition("/")[-1].rpartition(".NTD")[0]}.Clustered.NTD.fasta'

    cdhit_cmd  = 'cdhit-est -M 0 -G 0 -aL 0.0005 -aS 1.0 ' \
                f'-T {threads} ' \
                f'-c {pid} ' \
                f'-i {fasta_file} ' \
                f'-o {cdhit_out_fasta}'

    cdhit_rslt = subprocess.run(cdhit_cmd.split(),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True)

    return cdhit_out_fasta


def orf_calling(
        start_time,
        out_dir: str,
        taxon_name: str,
        fasta_file: str,
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
        blast_based: bool = True,
        top_hits_only: bool = False,
        prots: bool = False,
        verbose: bool = False) -> list:


    Path(out_dir).mkdir(exist_ok = True, parents = True)

    if blast_based:
        orf_call_fasta, og_hit_tsv = oc.orf_call_diamond(
                                            start_time,
                                            out_dir,
                                            fasta_file,
                                            taxon_name,
                                            diamond_db,
                                            delim,
                                            og_delim,
                                            gen_code,
                                            threads,
                                            evalue,
                                            min_id,
                                            max_hits,
                                            min_hit_cover,
                                            min_qry_cover,
                                            top_hits_only,
                                            prots,
                                            verbose)

        return orf_call_fasta, og_hit_tsv

    else:
        pass


def prep_transcriptomes(
        start_time,
        fasta_file: str,
        taxon_dir: str,
        taxon_name: str,
        taxon_code: str = None,
        orf_db: str = None,
        delim: str = '|',
        og_delim: str = '|',
        min_len: int = 200,
        gen_code: str = '1',
        evalue: float = 1e-10,
        min_id: int = 0,
        max_hits: int = 1,
        min_hit_cover: int = 60,
        min_qry_cover: int = 0,
        top_hits_only: bool = False,
        blast_based: bool = True,
        threads: int = 4,
        verbose: bool = False) -> None:

    out_dir_bu = f'{taxon_dir}/Original/'
    out_dir_ss = f'{taxon_dir}/Filter_Steps/Length_Filter/'
    out_dir_rr = f'{taxon_dir}/Filter_Steps/rRNA_Filter/'
    out_dir_oc = f'{taxon_dir}/OG_ORF_Calling/'
    out_dir_of = f'{taxon_dir}/Filter_Steps/Redundant_Transcript_Filter/'

    tax = taxon_name

    if taxon_code:
        tax = taxon_code

    if verbose:
        print('\n#------ Preparing Transcriptome Data -------#')
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Backing Up Data for {tax}')

    back_up_file(
        out_dir_bu,
        fasta_file)

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Removing short transcripts')

    sfilt_fasta = remove_short_seqs(
                    out_dir_ss,
                    taxon_name,
                    taxon_code,
                    fasta_file)

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Removing rRNA contamination')

    rrfilt_fasta = remove_rRNA(
                        out_dir_rr,
                        tax,
                        sfilt_fasta,
                        min_len,
                        threads)

    if verbose:
        print('\n#------ OG Assignment and ORF Calling ------#')

    init_og_fasta, og_hit_tsv = orf_calling(
                                    start_time,
                                    out_dir_oc,
                                    tax,
                                    rrfilt_fasta,
                                    orf_db,
                                    delim,
                                    og_delim,
                                    gen_code,
                                    threads,
                                    evalue,
                                    min_id,
                                    max_hits,
                                    min_hit_cover,
                                    min_qry_cover,
                                    blast_based,
                                    top_hits_only,
                                    verbose = verbose)

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Filtering Redundant ORFs')

    final_og_ntd_fasta = collapse_isoforms(
                            out_dir_of,
                            init_og_fasta,
                            threads = threads)

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Translating ORFs')

    final_og_aa_fasta = oc.translate_orfs(
                            final_og_ntd_fasta,
                            gen_code)

    if verbose:
        print('\n#--------- Backing Up Prepared ORFs --------#')

    back_up_file(
        taxon_dir,
        final_og_ntd_fasta)

    back_up_file(
        taxon_dir,
        final_og_aa_fasta)

    print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Back Up Finished')


def prep_wgs(
        start_time,
        fasta_file: str,
        taxon_dir: str,
        taxon_name: str,
        taxon_code: str = None,
        orf_db: str = None,
        delim: str = '|',
        og_delim: str = '|',
        gen_code: str = '1',
        evalue: float = 1e-10,
        min_id: int = 0,
        max_hits: int = 1,
        min_hit_cover: int = 60,
        min_qry_cover: int = 0,
        top_hits_only: bool = False,
        blast_based: bool = True,
        genbank: bool = True,
        prokaryotic: bool = False,
        remove_isoforms: bool = True,
        threads: int = 4,
        verbose: bool = False) -> None:


    out_dir_bu = f'{taxon_dir}/Original/'
    out_dir_oc = f'{taxon_dir}/OG_Calling/'

    tax = taxon_name

    if taxon_code:
        tax = taxon_code

    if verbose:
        print('\n#----- Preparing Annotated Genome Data -----#')
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Backing Up Data for {tax}')


    back_up_file(
        out_dir_bu,
        fasta_file)

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Preparing CDSs')

    if genbank:
        if verbose and remove_isoforms:
            print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Extracting Longest Isoforms')

        orf_fasta = eval_wgs_genbank(
                        start_time,
                        out_dir_bu,
                        taxon_name,
                        fasta_file,
                        taxon_code,
                        prokaryotic,
                        remove_isoforms,
                        delim,
                        verbose)

    else:
        orf_fasta = eval_wgs_non_gbk(
                        out_dir_bu,
                        taxon_name,
                        fasta_file,
                        taxon_code,
                        prokaryotic,
                        delim)

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Translating CDSs')

    peptide_fasta = oc.translate_orfs(
                    orf_fasta,
                    gen_code)

    if verbose:
        print('\n#------ OG Assignment and ORF Calling ------#')

    og_fasta, og_hit_tsv = orf_calling(
                                start_time,
                                out_dir_oc,
                                tax,
                                peptide_fasta,
                                orf_db,
                                delim,
                                og_delim,
                                gen_code,
                                threads,
                                evalue,
                                min_id,
                                max_hits,
                                min_hit_cover,
                                min_qry_cover,
                                blast_based,
                                top_hits_only,
                                True,
                                verbose)

    og_ntd_seqs = []
    og_pep_seqs = {i.id.rpartition(delim)[0]:i.id for i in SeqIO.parse(og_fasta,'fasta')}

    for i in SeqIO.parse(orf_fasta, 'fasta'):
        if i.id in og_pep_seqs:
            i.id = og_pep_seqs[i.id]
            i.description = ''
            i.name = ''
            og_ntd_seqs.append(i)

    SeqIO.write(og_ntd_seqs, og_fasta.replace(".AA.fasta",".NTD.fasta"),'fasta')

    if verbose:
        print('\n#-------- Backing Up Finalized Data --------#')

    back_up_file(
        taxon_dir,
        og_fasta)

    back_up_file(
        taxon_dir,
        og_fasta.replace(".AA.fasta",".NTD.fasta"))


def guess_data_type(
        start_time,
        fasta_file: str,
        taxon_dir: str,
        taxon_name: str,
        taxon_code: str = None,
        gen_code: str = '1',
        verbose = False) -> bool:
    """
    Guess the data-type based on frequency of in-frame stop codons.

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    taxon_dir:   Ouput direcotry to store data
    taxon_name:  taxon/sample name
    taxon_code:  abbreviated taxon/sample code
    gen_code:    ncbi translation table (determines stop codons)

    Returns
    ----------
    likely_transcripts:  whether the data are likely transcripts or ORFs
    """

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Evaluating data type ' \
            '(transcripts or CDSs) by in-frame stop codon frequency')

    out_dir_bu = f'{taxon_dir}/Original/'

    likely_transcripts = False

    tax = taxon_name

    if taxon_code:
        tax = taxon_code

    out_tsv = f'{out_dir_bu}{tax}.InFrame_Stop_Codons.tsv'

    back_up_file(
        out_dir_bu,
        fasta_file)

    stop_dict = {'1':['TGA','TAA','TAG', 'Total'], '4':['TAA','TAG', 'Total'],
        '6':['TGA', 'Total'], '10':['TAA','TAG', 'Total'], '25':['TAA','TAG', 'Total'],
        '29':['TGA', 'Total'], '30':['TGA', 'Total']}

    try:
        stop_codons = {i:0 for i in stop_dict[gen_code]}

    except KeyError:
        stop_codons = {i:0 for i in stop_dict['1']}

    stop_keys = list(stop_codons.keys())[:-1]

    for i in SeqIO.parse(fasta_file,'fasta'):
        for n in range(0, len(i) - 3, 3):
            codon = f'{i.seq[n:n+3]}'
            stop_codons['Total'] += 1
            if codon in stop_codons:
                stop_codons[codon] += 1

    for i in stop_keys:
        stop_freq = stop_codons[i]/(stop_codons['Total'] * 0.001)
        stop_codons[f'{i}_Norm'] = f'{stop_freq: .3f}'
        if float(stop_codons[f'{i}_Norm']) >= 1.5:
            likely_transcripts = True

    with open(out_tsv,'w+') as w:
        w.write('Stop_Codon\tCounts\tTotal_Codons\tFrequency_per_1k_Codons\n')
        for i in stop_keys:
            stop_norm = f'{i}_Norm'
            w.write(f'{i}\t{stop_codons[i]}\t{stop_codons["Total"]}\t{stop_codons[stop_norm]}\n')

    if verbose:
        if likely_transcripts:
            print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Data type: likely transcripts')

        else:
            print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Data type: likely CDSs')

    return likely_transcripts
