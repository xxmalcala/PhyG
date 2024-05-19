#!/usr/bin/env python3

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np

from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from ete3 import NCBITaxa

# Library Dependencies: BioPython, ete3
# Software Dependencies: CD-HIT, DIAMOND, EPA-ng, Gappa, MAFFT


def back_up_data(taxon_name: str, fasta_file: str) -> None:

    Path(f'{taxon_name}_PhyG_Contam/Original').mkdir(parents = True, exist_ok = True)
    shutil.copy2(fasta_file, f'{taxon_name}_PhyG_Contam/Original/')


def find_taxonomy(genus_species:str) -> list:
    ncbi = NCBITaxa()
    name2taxid = ncbi.get_name_translator([' '.join(genus_species.split('_')[:2])])

    if 'Spironema' in genus_species:
        name2taxid = ncbi.get_name_translator(['Spironema cf. multiciliatum'])
    elif 'vickermanii' in genus_species:
        name2taxid = ncbi.get_name_translator(['Mantamonas'])

    try:
        taxid = list(name2taxid.values())[0][0]
    except IndexError:
        taxid = None
    if taxid:
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        txnmy = [names[taxid] for taxid in lineage]
        return txnmy
    else:
        return None

def parse_taxonomy(genus_species: str) -> None:
    txnmy = find_taxonomy(genus_species)

    try:
        genus, species = genus_species.split('_')[:2]
        txn = genus[0] + species[:3]
    except ValueError:
        genus = genus_species.split('_')[0]
        species = 'sp'
        txn = genus[:2] + species

    if not txnmy:
        print('[Warning]: No NCBI taxonomy was found. Defaulting to generic taxon-code ' \
            f'"Un_un_{txn}"')
        return f'Un_un_{txn}'

    invalid_clades = set(['Bacteria','Archaea','Viruses'])

    if set(txnmy)&invalid_clades:
         print('[Warning]: Only supports Eukaryotes for now, sorry!')
         sys.exit(1)

    mjr_clades = {'Sar':'Sr','Discoba':'Ex','Metamonada':'Ex','Amoebozoa':'Am',
                'Viridiplantae': 'Pl', 'Rhodophyta': 'Pl', 'Opisthokonta': 'Op',
                'Glaucocystophyceae':'Pl'}

    mnr_clades = {'Sr':{'Ciliophora':'ci','Dinophyceae':'di','Apicomplexa':'ap',
                    'Stramenopila':'st','Rhizaria':'rh', 'Perkinsozoa':'pe'},
                'Ex':{'Euglenozoa':'eu','Heterolobosea':'he','Jakobida':'ja',
                    'Fornicata':'fo','Parabasalia':'pa','Preaxyostyla':'pr'},
                'Pl':{'Viridiplantae':'gr','Rhodophyta':'rh','Glaucocystophyceae':'gl'},
                'Op':{'Metazoa':'me','Choanoflagellata':'ch','Filasteria':'fi',
                    'Aphelida':'ap','Ichthyosporea':'ic','Rotosphaerida':'ro'},
                'Am':{'Discosea':'di','Archaemoebae':'ar','Eumycetozoa':'my','Tubulinea':'tu'},
                'EE':{'Haptophyta':'ha','Centroplasthelida':'ce','CRuMs':'cm',
                    'Ancyromonadida':'an','Breviatea':'br','Apusozoa':'ap',
                    'Cryptophyceae':'cr','Provora':'pr','Rhodelphea':'rh',
                    'Hemimastigophora':'he','Malawimonadida':'ma'}}

    mjr_c = ''
    mnr_c = ''
    for i in txnmy:
        if not mjr_c and i in mjr_clades.keys():
            mjr_c = mjr_clades[i]

    if not mjr_c:
        mjr_c = 'EE'

    for i in txnmy:
        if mjr_c in mnr_clades.keys() and not mnr_c:
            if i in mnr_clades[mjr_c].keys():
                mnr_c = mnr_clades[mjr_c][i]

    if mnr_c == '':
        mnr_c = 'is'

    return f'{mjr_c}_{mnr_c}_{txn}|{genus_species}'


def rename_seqs(taxon_name: str, taxon_code: str, fasta_file: str) -> str:
    trxp = 1
    renamed_dict = {}
    renamed_seqs = []

    outf = f'{taxon_name}_PhyG_Contam/Original/{taxon_name}'

    for i in SeqIO.parse(fasta_file,'fasta'):
        new_name = f'{taxon_code}|Transcript_{trxp}_Len_{len(i.seq)}'
        if 'cov' in i.id.lower():
            new_name += f'_Cov_{float(i.id.split("_cov_")[1].split("_")[0]):.2f}'
        renamed_dict[i.id] = new_name
        i.id = new_name
        i.description = ''
        i.name = ''
        renamed_seqs.append(i)
        trxp += 1

    with open(f'{outf}.Renamed_Seqs.tsv','w+') as w:
        w.write('Original_Name\tTemporary_Name\n')
        for k, v in renamed_dict.items():
            w.write(f'{k}\t{v}\n')

    SeqIO.write(renamed_seqs, f'{outf}.Query_Seqs.fasta','fasta')

    return f'{outf}.Query_Seqs.fasta'


def diamond_search(taxon_name: str, fasta_file: str, database: str, threads: int = 24, is_nuc: bool = True) -> str:
    # use a simple BLASTP command for now
    outf = f'{taxon_name}_PhyG_Contam/BLAST_Reports/{taxon_name}'
    dmnd_cmd = f'diamond blastp --quiet ' \
                f'-q {fasta_file} ' \
                f'-d {database} ' \
                f'-p {threads} ' \
                '-e 1e-10 ' \
                '-k 1 ' \
                '--min-orf 75 ' \
                '--subject-cover 60 ' \
                f'-o {outf}.BLASTP_Hits.tsv ' \
                '-f 6'
    if is_nuc:
        dmnd_cmd = dmnd_cmd.replace('blastp','blastx')
        dmnd_cmd += ' qseqid sseqid pident length qstart qend evalue bitscore qframe'

    Path(f'{taxon_name}_PhyG_Contam/BLAST_Reports').mkdir(parents = True, exist_ok = True)

    os.system(dmnd_cmd)

    return f'{outf}.BLASTP_Hits.tsv'


# IMPORTANT: NEED METHOD TO EXCISE THE ORF!!!

def assign_seq_ogs(taxon_name: str, fasta_file: str, og_hits_file: str, is_nuc: bool = True) -> str:
    clean_orfs = []
    updated_og_seqs = defaultdict(list)

    outf = f'{taxon_name}_PhyG_Contam/Query_ORFs/'

    tmp_seqs = {seq_rec.id:seq_rec for seq_rec in SeqIO.parse(fasta_file,'fasta')}

    for line in open(og_hits_file).readlines():
        txp, orf, qs, qe, qf = np.array(line.split('\t'))[[0,1,4,5,-1]]

        seq_rec = tmp_seqs[txp]

        seq_rec.id += f'{orf[-10:]}'
        seq_rec.description = ''
        seq_rec.name = ''

        if int(qf) > 0:
            seq_rec.seq = seq_rec.seq[int(qs)-1:int(qe)]

        elif int(qf) < 0:
            seq_rec.seq = seq_rec.seq[int(qe)-1:int(qs)].reverse_complement()

        updated_og_seqs[orf[-10:]].append(seq_rec)

    if is_nuc:
        nuc_outf = f'{taxon_name}_PhyG_Contam/Query_ORFs/Nucleotides/'
        Path(nuc_outf).mkdir(parents = True, exist_ok = True)

    for k, v in updated_og_seqs.items():
        SeqIO.write(v, f'{nuc_outf}{k}.Query_Seqs.NTD.fasta','fasta')

    if is_nuc:
        return nuc_outf

    return outf





# def assign_seq_ogs(taxon_name: str, aln_file: str, og_hits_file: str, is_nuc: bool = True) -> str:
#     # Note that this needs an out-put d
#     query_og_assignment = {}
#     updated_og_seqs = defaultdict(list)
#
#     outf = f'{taxon_name}_PhyG_Contam/Query_ORFs/'
#
#     for line in open(og_hits_file).readlines():
#         if is_nuc:
#             query_og_assignment[line.split('\t')[0]] = [line.split('\t')[1][-10:], int(line.split('\t')[-1])]
#         else:
#             query_og_assignment[line.split('\t')[0]] = [line.split('\t')[1][-10:]]
#
#     for seq_rec in SeqIO.parse(aln_file,'fasta'):
#         if seq_rec.id in query_og_assignment:
#             og = query_og_assignment[seq_rec.id][0]
#
#             if is_nuc:
#                 frame = query_og_assignment[seq_rec.id][-1]
#
#                 if frame > 0:
#                     tmp_seq = seq_rec.seq[frame-1:]
#
#                 elif frame < 0:
#                     tmp_seq = seq_rec.seq.reverse_complement()[abs(frame)-1:]
#
#                 end_trim = len(tmp_seq)%3
#                 seq_rec.seq = tmp_seq[:len(tmp_seq)-end_trim]
#
#             seq_rec.id += f'_{og}'
#             seq_rec.description = ''
#             seq_rec.name = ''
#
#             updated_og_seqs[og].append(seq_rec)
#
#     if is_nuc:
#         nuc_outf = f'{taxon_name}_PhyG_Contam/Query_ORFs/Nucleotides/'
#         Path(nuc_outf).mkdir(parents = True, exist_ok = True)
#
#     for k, v in updated_og_seqs.items():
#         SeqIO.write(v, f'{nuc_outf}{k}.Query_Seqs.NTD.fasta','fasta')
#
#     if is_nuc:
#         return nuc_outf
#
#    return outf


def collapse_redundant(query_fasta: str, threads: int = 24) -> str:
    outf = query_fasta.replace("Nucleotides","Clustered").replace("fasta","97idClust.fasta")

    cd_hit_cmd = f'cd-hit-est ' \
                '-G 0 ' \
                '-c 0.97 ' \
                '-aS 1.0 ' \
                '-aL .005 ' \
                '-M 4000 ' \
                f'-T {threads} ' \
                f'-i {query_fasta} ' \
                f'-o {outf}'

    cdht_rslt = subprocess.run(
                    cd_hit_cmd.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True
                    )

    return outf


def clust_nuc_fastas(taxon_name: str, nuc_dir: str, threads: int = 24) -> list:
    clust_files = []

    Path(f'{taxon_name}_PhyG_Contam/Query_ORFs/Clustered/').mkdir(parents = True, exist_ok = True)

    for f in glob.glob(f'{nuc_dir}*.NTD.fasta'):
        clust_files.append(collapse_redundant(f, threads))

    return clust_files


def translate_seqs(nuc_fasta_file: str, gcode: int = 1) -> None:
    pep_seqs = []
    for i in SeqIO.parse(nuc_fasta_file,'fasta'):
        i.seq = i.seq.translate(gcode).rstrip('*').replace("*","X")
        pep_seqs.append(i)

    outf = nuc_fasta_file.replace("Clustered/","Peptides/").rpartition(".NTD.")[0]

    SeqIO.write(pep_seqs, f'{outf}.AA.fasta','fasta')


def prep_for_epa_ng(taxon_code: str, pep_dir: str, ref_msa_dir: str, ref_tree_dir: str, threads: int = 24) -> None:
    msa_num = {}

    query_msa_dir = pep_dir.replace("Peptides","Query_MSAs")

    Path(query_msa_dir).mkdir(parents = True, exist_ok = True)

    for f in glob.glob(f'{pep_dir}*AA.fasta'):
        msa_num[f.rpartition('/')[-1].partition('.')[0]] = [f]

    for f in glob.glob(f'{ref_msa_dir}/*fasta'):
        og_num = f.rpartition('/')[-1].partition('.')[0]
        if og_num in msa_num.keys():
            msa_num[og_num].append(f)

    for f in glob.glob(f'{ref_tree_dir}/*'):
        og_num = f.rpartition('/')[-1].partition('.')[0]
        if og_num in msa_num.keys():
            msa_num[og_num].append(f)

    for k, v in msa_num.items():
        query_msa = mafft_align(taxon_code, query_msa_dir, v[1], v[0], threads)
        msa_num[k].append(query_msa)

    return msa_num


def mafft_align(taxon_code: str, outdir: str, ref_msa: str, query_seqs: str, threads: int = 24):
    # dirty for now...
    outf = f'{outdir}{query_seqs.rpartition("/")[-1].replace("fasta","MAFFT.fasta")}'

    mafft_cmd = f'mafft --quiet ' \
                f'--thread {threads} ' \
                '--auto ' \
                f'--add {query_seqs} ' \
                f'--keeplength {ref_msa} ' \
                f'> {outf}'


    os.system(mafft_cmd)

    # DEBUG: subprocess not working for some reason...

    # mafft_rslt = subprocess.run(
    #                 mafft_cmd.split(),
    #                 stdout=subprocess.PIPE,
    #                 stderr=subprocess.PIPE,
    #                 universal_newlines=True
    #                 )

    query_msa = [i for i in SeqIO.parse(outf,'fasta') if taxon_code in i.id]
    SeqIO.write(query_msa, outf, 'fasta')

    return outf


def run_epa_ng(taxon_name: str, epa_dict: dict, threads: int = 24):
    outdir = f'{taxon_name}_PhyG_Contam/Updated_Trees/'

    Path(outdir).mkdir(parents = True, exist_ok = True)

    for k, v in epa_dict.items():
        if len(v) != 4:
            print(k, v)
            sys.exit()
        ref_msa = v[1]
        ref_tree = v[2]
        query_msa = v[3]

        epa_ng_cmd = f'epa-ng --redo ' \
                        f'-T {threads} ' \
                        f'--ref-msa {ref_msa} ' \
                        f'--tree {ref_tree} ' \
                        f'--query {query_msa} ' \
                        f'-w {outdir} ' \
                        '--model LG+I+G4'

        epang_rslt = subprocess.run(
                        epa_ng_cmd.split(),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True
                        )
        out_jplace = f'{outdir}{k}.Updated.jplace'
        shutil.copy2(f'{outdir}epa_result.jplace',f'{out_jplace}')

        convert_jplace_newick(outdir, out_jplace)

        os.system(f'rm {outdir}epa_*')
        os.system(f'rm {outdir}*jplace')


def convert_jplace_newick(outdir: str, jplace_file: str):
    gappa_cmd = f'gappa examine graft ' \
                '--fully-resolve ' \
                f'--jplace-path {jplace_file} ' \
                f'--out-dir {outdir}'

    gappa_rslt = subprocess.run(
                    gappa_cmd.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True
                    )


if __name__ == '__main__':
    try:
        fasta_file = sys.argv[1]
        database = sys.argv[2]
        taxon_name = sys.argv[3]
        ref_msa_dir = sys.argv[4]
        ref_tree_dir = sys.argv[5]

    except:
        print('Usage:\n\n    python3 phyg_contam.py [FASTA-FILE] [DATABASE] [TAXON-NAME] '
                '[REF-MSA-DIR] [REF-TREE-DIR]\n')
        sys.exit(1)

    ## check that valid taxon-code is NOT given as the taxon-name...

    back_up_data(taxon_name, fasta_file)

    taxon_code = parse_taxonomy(taxon_name)

    query_fasta = rename_seqs(taxon_name, taxon_code, fasta_file)

    og_hits_file = diamond_search(taxon_name, query_fasta, database)

    # IMPORTANT: NEED METHOD TO EXCISE THE ORF!!!

    dir_to_clust = assign_seq_ogs(taxon_name, query_fasta, og_hits_file)

    nuc_to_pep = clust_nuc_fastas(taxon_name, dir_to_clust)

    pep_dir = f'{taxon_name}_PhyG_Contam/Query_ORFs/Peptides/'

    Path(pep_dir).mkdir(parents = True, exist_ok = True)

    for nuc_fasta in nuc_to_pep:
        translate_seqs(nuc_fasta)

    data_for_epa = prep_for_epa_ng(taxon_code, pep_dir, ref_msa_dir, ref_tree_dir)
    run_epa_ng(taxon_name, data_for_epa)
