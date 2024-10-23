#!/usr/bin/env python3

import logging, shutil, subprocess, sys

from collections import defaultdict

from pathlib import Path

import largest_isoform as liso
import phyg_orf_related as oc

from Bio import SeqIO

"""
notes to self:

if given ORFs from WGS file (from GenBank/RefSeq), must extract the accession codes and can
ignore the minimum length requirements

Similarly, if given a file of ORFs, no need to manipulate the sequences inside...
-- just add the taxon_name and delimiter IF the seq-id doesn't start with taxon_name
"""

def back_up_file(
        out_dir: str,
        fasta_file: str) -> None:

    Path(out_dir).mkdir(exist_ok = True, parents = True)

    shutil.copy2(fasta_file, out_dir)


def eval_wgs_genbank(
        out_dir: str,
        taxon_name: str,
        fasta_file: str,
        taxon_code: str = None,
        prokaryotic: bool = False,
        remove_isoforms: bool = True,
        delim: str = '|') -> str:

    Path(out_dir).mkdir(exist_ok = True, parents = True)


    final_filt_fas = f'{out_dir}{taxon_name}.LongestIsoform'

    if taxon_code:
        final_filt_fas = final_filt_fas.replace(taxon_name, taxon_code)

    if not remove_isoforms:
        final_filt_fas.replace(".LongestIsoform",".ORFs")

    final_seqs = []
    name_conversion = {}

    if remove_isoforms:
        isoforms_to_keep = liso.save_isoforms(fasta_file, True, prokaryotic)
    else:
        isoforms_to_keep = liso.save_isoforms(fasta_file, False, prokaryotic)

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
        delim: str = '|'
        ) -> str:


    Path(out_dir).mkdir(exist_ok = True, parents = True)

    final_filt_fas = f'{out_dir}{taxon_name}.ORFs.fasta'

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

    SeqIO.write(final_seqs, final_filt_fas, 'fasta')

    with open(final_filt_fas.replace(".fasta",".NameConversion.tsv"), 'w+') as w:
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
    transcript:

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
    taxon_name:  species/taxon name or abbreviated code
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
            prots: bool = False) -> str:


    Path(out_dir).mkdir(exist_ok = True, parents = True)

    if blast_based:
        orf_call_fasta, og_hit_tsv = oc.orf_call_diamond(
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
                                            prots)

        return orf_call_fasta, og_hit_tsv
    else:
        pass


def prep_transcriptomes(
        taxon_dir: str,
        taxon_name: str,
        taxon_code: str = None,
        orf_db: str = None,
        delim: str = '|',
        og_delim: str = 'OG6',
        min_len: int = 200,
        gen_code: str = '1',
        evalue: float = 1e-10,
        min_id: int = 0,
        max_hits: int = 1,
        min_hit_cover: int = 60,
        min_qry_cover: int = 0,
        top_hits_only: bool = False,
        blast_based: bool = True,
        threads: int = 4):

    out_dir_bu = f'{taxon_dir}/Original/'
    out_dir_ss = f'{taxon_dir}/Filter_Steps/Length_Filter/'
    out_dir_rr = f'{taxon_dir}/Filter_Steps/rRNA_Filter/'
    out_dir_oc = f'{taxon_dir}/OG_ORF_Calling/'
    out_dir_of = f'{taxon_dir}/Filter_Steps/Redundant_Transcript_Filter/'

    tax = taxon_name

    if taxon_code:
        tax = taxon_code

    back_up_file(
        out_dir_bu,
        fasta_file
        )

    sfilt_fasta = remove_short_seqs(
                    out_dir_ss,
                    taxon_name,
                    taxon_code,
                    fasta_file
                    )

    rrfilt_fasta = remove_rRNA(
                        out_dir_rr,
                        tax,
                        sfilt_fasta,
                        min_len,
                        threads
                        )

    init_og_fasta, og_hit_tsv = orf_calling(
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
                                    top_hits_only
                                    )

    final_og_ntd_fasta = collapse_isoforms(
                            out_dir_of,
                            init_og_fasta,
                            threads = threads
                            )

    final_og_aa_fasta = oc.translate_orfs(
                            final_og_ntd_fasta,
                            gen_code
                            )

    back_up_file(
        taxon_dir,
        final_og_ntd_fasta
        )

    back_up_file(
        taxon_dir,
        final_og_aa_fasta
        )


def prep_wgs(
        taxon_dir: str,
        taxon_name: str,
        taxon_code: str = None,
        orf_db: str = None,
        delim: str = '|',
        og_delim: str = 'OG6',
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
        threads: int = 4):

    out_dir_bu = f'{taxon_dir}/Original/'
    out_dir_oc = f'{taxon_dir}/OG_Calling/'

    tax = taxon_name

    if taxon_code:
        tax = taxon_code

    back_up_file(
        out_dir_bu,
        fasta_file
        )

    genbank = True

    if genbank:
        orf_fasta = eval_wgs_genbank(
                        out_dir_bu,
                        taxon_name,
                        fasta_file,
                        taxon_code,
                        prokaryotic,
                        remove_isoforms,
                        delim
                        )

    else:
        orf_fasta = eval_wgs_non_gbk(
                        out_dir_bu,
                        taxon_name,
                        fasta_file,
                        taxon_code,
                        prokaryotic,
                        delim
                        )

    peptide_fasta = oc.translate_orfs(
                    orf_fasta,
                    gen_code
                    )

    og_fasta, og_hit_tsv = orf_calling(
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
                                True
                                )

    og_ntd_seqs = []
    og_pep_seqs = {i.id.rpartition(delim)[0]:i.id for i in SeqIO.parse(og_fasta,'fasta')}
    for i in SeqIO.parse(orf_fasta, 'fasta'):
        if i.id in og_pep_seqs:
            i.id = og_pep_seqs[i.id]
            i.description = ''
            i.name = ''
            og_ntd_seqs.append(i)

    SeqIO.write(og_ntd_seqs, og_fasta.replace(".AA.fasta",".NTD.fasta"),'fasta')

    back_up_file(
        taxon_dir,
        og_fasta
        )

    back_up_file(
        taxon_dir,
        og_fasta.replace(".AA.fasta",".NTD.fasta")
        )

if __name__ == '__main__':
    try:
        fasta_file = sys.argv[1]
        taxon_name = sys.argv[2]
        orf_db = sys.argv[3]
    except:
        print('Usage:\n    python3 phyg_add_taxa.py [FASTA-FILE] [TAXON-NAME/CODE] [ORF-DATABASE]')
        sys.exit(1)

    taxon_code = None # 'Sr_ci_Samb'

    taxon_dir = f'{taxon_name}_PhyG'

    if taxon_code:
        taxon_dir = taxon_dir.replace(taxon_name, taxon_code)

    transcriptome = False
    wgs = True

    if transcriptome:
        prep_transcriptomes(
                taxon_dir,
                taxon_name,
                taxon_code,
                orf_db,
                threads = 26
                )

    elif wgs:
        genbank = True
        remove_isoforms = True
        prokaryotic = True
        prep_wgs(
                taxon_dir,
                taxon_name,
                taxon_code,
                orf_db,
                genbank = genbank,
                prokaryotic = prokaryotic,
                remove_isoforms = remove_isoforms,
                threads = 26
                )
