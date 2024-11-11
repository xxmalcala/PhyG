#!/usr/bin/env python3

import logging, shutil, subprocess, sys

from collections import defaultdict

from pathlib import Path

# from bin import largest_isoform as liso
from bin import phyg_orf_related as oc

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


"""
notes to self:

if given ORFs from WGS file (from GenBank/RefSeq), must extract the accession codes and can
ignore the minimum length requirements

Similarly, if given a file of ORFs, no need to manipulate the sequences inside...
-- just add the taxon_name and delimiter IF the seq-id doesn't start with taxon_name
"""

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
        fasta_file: str,
        prok: bool = False) -> dict:
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


def save_isoforms(
        fasta_file: str,
        longest_isoform: bool = True,
        prok: bool = False) -> list:
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
        isoforms_to_keep = save_isoforms(fasta_file, True, prokaryotic)
        if not isoforms_to_keep:
            print('ERROR: Unable to identify NCBI loci and remove redundant isoforms'
            f' in {fasta_file.rpartition("/")[-1]}. Please ensure that this file was sourced'
            ' from GenBank/RefSeq.')
            print("If this warning persists, please open an issue on PhyG's GitHub page.")
            sys.exit(1)

    else:
        isoforms_to_keep = save_isoforms(fasta_file, False, prokaryotic)

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
        fasta_file: str,
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
        fasta_file: str,
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

    transcriptome = True
    wgs = False

    codes = {'Furgasonia':'6','Nassula': '6', 'Mesodinium': '29'}
    gcode = '1'

    if taxon_name.split("_")[0] in codes:
        gcode = codes[taxon_name.split("_")[0]]

    if transcriptome:
        prep_transcriptomes(
                fasta_file,
                taxon_dir,
                taxon_name,
                taxon_code,
                orf_db,
                og_delim = '|',
                gen_code = gcode,
                threads = 26
                )

    elif wgs:
        genbank = True
        remove_isoforms = True
        prokaryotic = True
        prep_wgs(
            fasta_file,
            taxon_dir,
            taxon_name,
            taxon_code,
            orf_db,
            genbank = genbank,
            prokaryotic = prokaryotic,
            remove_isoforms = remove_isoforms,
            threads = 26
            )
