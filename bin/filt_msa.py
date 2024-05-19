#!/usr/bin/env python3

"""
Controls the flow for Guidance2 generation, assessment, and filtering.
"""

import glob, logging, os, shutil, subprocess, sys, time
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


def pre_filt(fasta_file):
    seqs = [len(i.seq) for i in SeqIO.parse(fasta_file,'fasta')]
    return sum(seqs)/len(seqs)


def clust_seqs(fasta_file, seq_id = 95, threads = 4):
    if seq_id > 1:
        seq_id = seq_id/100

    clust_cmd = f'cd-hit ' \
                '-G 0 ' \
                '-aS 0.67 ' \
                '-aL 0.005 ' \
                '-g 1 ' \
                '-d 0 ' \
                f'-T {threads} ' \
                f'-c {seq_id} ' \
                f'-i {fasta_file} ' \
                f'-o {fasta_file.replace(".Guid1",".clust.Guid1")}'

    clust_call = subprocess.call(clust_cmd, shell = True,
        stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

    return f'{fasta_file.replace(".Guid1",".clust.Guid1")}.clstr'


def prep_for_guid(fasta_file, guid_dir, guidance_iter = 1):
    # Create temporary fasta files to track the Guidance2 Iterations
    pg_fas = f'{fasta_file.split("/")[-1].split(".fa")[0]}.Guid{guidance_iter}.fas'
    seq_codes = {}
    num_seqs = []
    seq_count = 1

    if guidance_iter == 1:
        og_mean_len = pre_filt(fasta_file)

    # Rename/number the sequences for Guidance2
    for i in SeqIO.parse(fasta_file,'fasta'):
        if guidance_iter == 1:
            if 0.5*og_mean_len <= len(i.seq) <= 2*og_mean_len:
                num_seqs.append(f'>{seq_count}\n{i.seq.replace("*","X")}')
                seq_codes[f'{seq_count}'] = i.description
                seq_count += 1
        else:
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

    if guidance_iter == 1:
        # print('here')
        # clust_file = clust_seqs(f'{guid_dir}{pg_fas}', seq_id = 95, threads = 4)
        # print(clust_file.rpartition(".")[0], f'{guid_dir}{pg_fas}')
        return f'{guid_dir}{pg_fas}', f'{guid_dir}{pg_fas}', seq_codes

    return f'{guid_dir}{pg_fas}', seq_codes


def run_guidance(guidance_path: str,
                tmp_msa: str,
                seq_cutoff: float = 0.6,
                col_cutoff: float = 0.93,
                threads: int = 1,
                iteration: int = 1,
                quiet: bool = False) -> str:

    temp_outdir = tmp_msa.rstrip('.fas')
    # print(tmp_msa)

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
    if open(msa_file).read().count('>') > 3:
        return True

    else:
        return False


def backup_guid_files(tmp_outdir, scdir, seq_codes):
    pfx = tmp_outdir.split('/')[-1]
    orig_ssc_file = f'{tmp_outdir}/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names'
    orig_csc_file = f'{tmp_outdir}/MSA.MAFFT.Guidance2_res_pair_col.scr.csv'
    html_file = f'{tmp_outdir}/MSA.MAFFT.Guidance_res_pair_res.html'

    if not Path(orig_ssc_file).is_file():
        return None, None

    with open(f'{scdir}{pfx}.SeqScores.tsv', 'w+') as fout:
        fout.write('SeqName\tSeq_Score\n')
        for line in open(orig_ssc_file).readlines()[1:-1]:
            n = seq_codes[line.split('\t')[0]]
            s = line.split("\t")[1]
            fout.write(f'{n}\t{s}')

    # shutil.copy2(orig_ssc_file, f'{scdir}{pfx}.SeqScores.tsv')
    # shutil.copy2(clust_file, f'{scdir}{pfx}.Cluster')
    shutil.copy2(orig_csc_file, f'{scdir}{pfx}.ColScores.csv')
    shutil.copy2(html_file, f'{scdir}{pfx}.Scores.html')

    return orig_ssc_file, orig_csc_file

def eval_seqs(tmp_msa, seq_scores, seq_cutoff) -> list:
    return [i.split('\t')[0] for i in open(seq_scores).readlines()[1:-1] if float(i.split('\t')[1]) >= float(seq_cutoff)]


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


def parse_clust_file(seq_codes, clust_file):
    clust_summary = {}
    for cluster in open(clust_file).read().split('>Cluster '):
        if '*' in cluster:
            cl_seqs = []
            ref_seq = None
            for i in cluster.split('>'):
                if '*' in i:
                    ref_seq = i.split("...")[0]
                elif '...' in i:
                    cl_seqs.append(i.split("...")[0])
            if ref_seq:
                clust_summary[ref_seq] = cl_seqs

    return clust_summary


def add_seqs_to_msa(tmp_dir, msa_with_cols, seqs_to_add, fasta_file):

    ref_msa = f'{tmp_dir}/ref_msa.fas'
    qry_msa = f'{tmp_dir}/missing_msa.fas'
    fin_msa = f'{tmp_dir}/prepped_msa.fas'

    with open(f'{tmp_dir}/ref_msa.fas','w+') as w:
        for i in msa_with_cols:
            w.write(f'>{i[0]}\n{i[1]}\n')

    with open(f'{tmp_dir}/missing_msa.fas','w+') as w:
        for i in SeqIO.parse(fasta_file, 'fasta'):
            if i.id in seqs_to_add:
                w.write(f'>{i.id}\n{i.seq}\n')

    mafft_cmd = f'mafft --keeplength --add {qry_msa} {ref_msa} > {fin_msa}'

    mafft_call = subprocess.call(mafft_cmd, shell = True,
                stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

    final_msa = [(i.description, f'{i.seq}') for i in SeqIO.parse(fin_msa, 'fasta')]

    return final_msa


def finalize_guid_output(tmp_outdir, seq_codes, fasta_file, check_clust = False):
    fin_dir = f'{tmp_outdir.split("Guidance_")[0]}/Final_MSAs/'
    Path(fin_dir).mkdir(parents = True, exist_ok = True)

    fin_msa = f'{tmp_outdir.split("/")[-1].split(".Guid")[0]}.PostGuid.fas'
    fin_msa_cols_rm = f'{fin_msa.replace(".fas",".ColsRemoved.fas")}'

    msa_with_cols = [(i.description, f'{i.seq}') for i in SeqIO.parse(f'{tmp_outdir}/MSA.MAFFT.aln.With_Names', 'fasta')]

    msa_seqs = [i[0] for i in msa_with_cols]

    if check_clust:
        clust_summary = parse_clust_file(seq_codes, clust_file)
        seqs_to_add = []

        for k, v in clust_summary.items():
            if k in msa_seqs:
                seqs_to_add += v

        if seqs_to_add:

            fin_msa_with_cols = add_seqs_to_msa(tmp_outdir.rpartition("/")[0], msa_with_cols, seqs_to_add, fasta_file)

            rename_seqs(fin_msa_with_cols, seq_codes, f'{fin_dir}{fin_msa}')

        else:
            rename_seqs(msa_with_cols, seq_codes, f'{fin_dir}{fin_msa}')

    else:
        rename_seqs(msa_with_cols, seq_codes, f'{fin_dir}{fin_msa}')
    # rename_seqs(msa_cols_rm, seq_codes, f'{fin_dir}{fin_msa_cols_rm}')


def guid_eval(tmp_msa, seq_codes, tmp_outdir, scdir, fdir, seq_cutoff, max_iter, cur_iter):
    score_files = backup_guid_files(tmp_outdir, scdir, seq_codes)
    if not score_files:
        return None

    seqs_passing = eval_seqs(tmp_msa, score_files[0], seq_cutoff)

    if cur_iter <= max_iter:
        filt_fas = f'{fdir}{tmp_msa.split("/")[-1].replace(".fas",".RemovedSeqs.fas")}'
        pass_fas = f'{tmp_msa.split(".Guid")[0]}.Guid{cur_iter+1}.fas'

        if save_seqs(seqs_passing, tmp_msa, tmp_outdir, pass_fas, filt_fas):
            return pass_fas, True

        else:
            return tmp_msa, False


def iter_guid(guid_path: str,
                fasta_file: str,
                gdir: str,
                scdir: str,
                fdir: str,
                max_iter: int,
                threads: int,
                seq_cutoff: float,
                col_cutoff: float,
                quiet: bool = False) -> str:

    check_clust = False
    sttime = time.time()

    cur_tmp_fas, pg_fas, snames = prep_for_guid(fasta_file, gdir, 1)

    if max_iter == 0:
        pass

    if not check_fasta(cur_tmp_fas):
        check_clust = False
        cur_tmp_fas = pg_fas

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

        cur_tmp_fas, iter_again = guid_eval(cur_tmp_fas,
                                    snames,
                                    tmp_outdir,
                                    scdir,
                                    fdir,
                                    seq_cutoff,
                                    max_iter,
                                    n + 1)

        if not iter_again or n + 1 == max_iter:

            finalize_guid_output(tmp_outdir, snames, pg_fas, check_clust)

            if not quiet:
                print(f'     [{timedelta(seconds=round(time.time()-sttime))}] Guidance iterations finished\n')
            break


def clean_up(gdir):
    subprocess.call(f'rm -r {gdir}*.Guid*', shell = True,
                stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

    subprocess.call(f'rm -r {gdir}*.fas*', shell = True,
                stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

# os.system(f'rm -r {gdir}/*.Guid*')
    # os.system(f'rm {gdir}/*.fas*')


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

    init_fastas.sort(key = os.path.getsize)
    # msa_sizes = [(f, open(f).read().count('>')) for f in init_fastas]
    # msa_sizes.sort(key=lambda x: x[1])
    #
    # init_fastas = [i[0] for i in msa_sizes]


    if not quiet:
        print('#--- MSA Building and Homology Assessment ---#')

    for fas in init_fastas:
        if not quiet:
                    #--- MSA Building and Homology Assessment ---#
            print(f'|--- {fas.split("/")[1].split(".fa")[0]}')

        iter_guid(guid_path,
                    fas,
                    gdir,
                    scdir,
                    fdir,
                    max_iter,
                    threads,
                    seq_cutoff,
                    col_cutoff,
                    quiet)

        if not keep_all:
            clean_up(gdir)

    return f'{gdir}/Final_MSAs/'
