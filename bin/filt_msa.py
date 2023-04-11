#!/usr/bin/env python3

#__author__: Xyrus X. Maurer-Alcala
#__last_updated__: 2023-04-10

""" Controls the flow for Guidance2 generation, assessment, and filtering."""

import glob, logging, shutil, subprocess, sys, time
import numpy as np
from Bio import SeqIO
from pathlib import Path
from datetime import timedelta



def prep_dir(proj_name):
    cwd = f'{Path.cwd()}/'
    gdir = f'{cwd}{proj_name}_PhyG_Results/Guidance_Steps/'
    sdir = f'{cwd}{proj_name}_PhyG_Results/Guidance_Steps/Guidance_Scores/'
    fdir = sdir.replace("e_Scores/","e_Filtered_Seqs/")

    Path(gdir).mkdir(parents = True, exist_ok = True)
    Path(sdir).mkdir(parents = True, exist_ok = True)
    Path(fdir).mkdir(parents = True, exist_ok = True)
    Path(f'{gdir}Converted_SeqCodes/').mkdir(parents = True, exist_ok = True)
    return gdir, sdir, fdir


def prep_for_guid(fasta_file, guid_dir, guidance_iter = 1):
    # Create temporary fasta files to track the Guidance2 Iterations
    pg_fas = f'{fasta_file.split("/")[-1].split(".fa")[0]}.Guid{guidance_iter}.fas'
    seq_codes = {}
    num_seqs = []
    seq_count = 1

    # Rename/number the sequences for Guidance2
    for i in SeqIO.parse(fasta_file,'fasta'):
        num_seqs.append(f'>{seq_count}\n{i.seq.replace("*","X")}')
        seq_codes[f'{seq_count}'] = i.description
        seq_count += 1

    # Write out the temporary files for the current iteration of Guidance2.
    # The aim is to protect the original FASTA files.
    with open(f'{guid_dir}{pg_fas}','w+') as w:
        w.write('\n'.join(num_seqs))

    with open(f'{guid_dir}Converted_SeqCodes/{pg_fas.split(".Guid")[0]}.SeqCodes.tsv','w+') as w:
        w.write('Original_Name\tTemporary_Number\n')
        for k, v in seq_codes.items():
            w.write(f'{v}\t{k}\n')

    return f'{guid_dir}{pg_fas}', seq_codes


def run_guidance(guidance_path: str,
                tmp_msa: str,
                seq_cutoff: float = 0.6,
                col_cutoff: float = 0.93,
                threads: int = 1,
                iteration: int = 1,
                quiet: bool = False) -> str:
    temp_outdir = tmp_msa.rstrip('.fas')
    guid_cmd = (f'perl {guidance_path} ' \
        f'--seqFile {tmp_msa} ' \
        f'--msaProgram MAFFT ' \
        f'--seqType aa ' \
        f'--outDir {temp_outdir} ' \
        f'--seqCutoff {seq_cutoff} ' \
        f'--colCutoff {col_cutoff} '
        f'--bootstraps 10 ' \
        f'--proc_num {threads} ' \
        f'--outOrder as_input ' \
        f'--MSA_Param \\\-"-auto ' \
        f'--maxiterate 200 '
        f'--thread {threads} '\
        f'--bl 62 --anysymbol"')
    guid_call = subprocess.call(guid_cmd, shell = True,
        stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
    return temp_outdir


def check_fasta(msa_file):
    if open(msa_file).read().count('>') > 4:
        return True
    else:
        return False


def backup_guid_files(tmp_outdir, scdir):
    pfx = tmp_outdir.split('/')[-1]
    orig_ssc_file = f'{tmp_outdir}/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names'
    orig_csc_file = f'{tmp_outdir}/MSA.MAFFT.Guidance2_res_pair_col.scr.csv'
    html_file = f'{tmp_outdir}/MSA.MAFFT.Guidance_res_pair_res.html'

    shutil.copy2(orig_ssc_file, f'{scdir}{pfx}.SeqScores.tsv')
    shutil.copy2(orig_csc_file, f'{scdir}{pfx}.ColScores.csv')
    shutil.copy2(html_file, f'{scdir}{pfx}.Scores.html')

    return orig_ssc_file, orig_csc_file


def eval_seqs(tmp_msa, seq_scores, seq_cutoff) -> list:
    return [i.split('\t')[0] for i in open(seq_scores).readlines()[1:-1] if float(i.split('\t')[1]) >= seq_cutoff]


def save_seqs(seq_list, tmp_msa, tmp_outdir, pass_seq_fas, filt_seq_fas):
    seqs_passing, seqs_removed = [], []

    for i in SeqIO.parse(tmp_msa,'fasta'):
        if i.description in seq_list:
            seqs_passing.append(f'>{i.description}\n{i.seq}\n')

        else:
            seqs_removed.append(f'>{i.description}\n{i.seq}\n')

    if seqs_removed:
        with open(f'{pass_seq_fas}','w+') as w:
            w.write(''.join(seqs_passing))

        with open(f'{filt_seq_fas}','w+') as w:
            w.write(''.join(seqs_removed))

        return True

    else:
        return False


def rename_seqs(seq_list, seq_codes, outfas):
    renamed_seqs = [f'>{seq_codes[i[0]]}\n{i[1]}\n' for i in seq_list]
    with open(outfas,'w+') as w:
        w.write(''.join(renamed_seqs))


def finalize_guid_output(tmp_outdir, seq_codes):
    fin_dir = f'{tmp_outdir.split("Guidance_")[0]}/Final_MSAs/'
    Path(fin_dir).mkdir(parents = True, exist_ok = True)

    fin_msa = f'{tmp_outdir.split("/")[-1].split(".Guid")[0]}.PostGuid.fas'
    fin_msa_cols_rm = f'{fin_msa.replace(".fas",".ColsRemoved.fas")}'

    msa_with_cols = [(i.description, f'{i.seq}') for i in SeqIO.parse(f'{tmp_outdir}/MSA.MAFFT.aln.With_Names', 'fasta')]
    msa_cols_rm = [(i.description, f'{i.seq}') for i in SeqIO.parse(f'{tmp_outdir}/MSA.MAFFT.Without_low_SP_Col.With_Names', 'fasta')]

    rename_seqs(msa_with_cols, seq_codes, f'{fin_dir}{fin_msa}')
    rename_seqs(msa_cols_rm, seq_codes, f'{fin_dir}{fin_msa_cols_rm}')


def guid_eval(tmp_msa, tmp_outdir, scdir, fdir, seq_cutoff, max_iter, cur_iter):
    score_files = backup_guid_files(tmp_outdir, scdir)
    seqs_passing = eval_seqs(tmp_msa, score_files[0], seq_cutoff)

    if cur_iter <= max_iter:
        filt_fas = f'{fdir}{tmp_msa.split("/")[-1].replace(".fas",".RemovedSeqs.fas")}'
        pass_fas = f'{tmp_msa.split(".Guid")[0]}.Guid{cur_iter+1}.fas'

        if save_seqs(seqs_passing, tmp_msa, tmp_outdir, pass_fas, filt_fas):
            return pass_fas
        else:
            return None


def iter_guid(guid_path: str,
                fasta_file: str,
                gdir: str,
                max_iter: int,
                threads: int,
                seq_cutoff: float,
                col_cutoff: float,
                quiet: bool = False) -> str:

    tmp_fas, snames = prep_for_guid(fasta_file, gdir, 1)

    if max_iter == 0:
        pass

    sttime = time.time()
    cur_tmp_fas = tmp_fas

    for n in range(max_iter):
        if not check_fasta(cur_tmp_fas):
            break
        if not quiet:
            print(f'     [{timedelta(seconds=round(time.time()-sttime))}] Guidance2 iteration [{n + 1}]')

        tmp_outdir = run_guidance(guid_path,
                            cur_tmp_fas,
                            seq_cutoff,
                            col_cutoff,
                            threads,
                            n+1,
                            quiet)

        cur_tmp_fas = guid_eval(cur_tmp_fas,
                            tmp_outdir,
                            scdir,
                            fdir,
                            seq_cutoff,
                            max_iter,
                            n + 1)

        if not cur_tmp_fas or n + 1 == max_iter:
            print('here')
            finalize_guid_output(tmp_outdir, snames)
            if not quiet:
                print(f'     [{timedelta(seconds=round(time.time()-sttime))}] Guidance iterations finished\n')


def clean_up(gdir):
    os.system(f'rm -r {gdir}/*.Guid*')


def make_msa_files(guid_path,
                    msa_dir: str,
                    proj_name: str,
                    max_iter: str = 10,
                    threads: str = 1,
                    seq_cutoff: float = 0.6,
                    col_cutoff: float = 0.7,
                    keep_all: bool = False,
                    quiet: bool = False):

    gdir, scdir, fdir = prep_dir(proj_name)


    init_fastas = glob.glob(f'{msa_dir}/*.fa*')

    if not quiet:
        print('#--- MSA Building and Homology Assessment ---#')

    for fas in init_fastas:
        if not quiet:
                    #--- MSA Building and Homology Assessment ---#
            print(f'|--- {fas.split("/")[1].split(".fa")[0]}')
        iter_guid(guid_path,
                    fas,
                    gdir,
                    max_iter,
                    threads,
                    seq_cutoff,
                    col_cutoff,
                    quiet)
        if not keep_all:
            clean_up(gdir)

    return f'{gdir}/Final_MSAs/'
