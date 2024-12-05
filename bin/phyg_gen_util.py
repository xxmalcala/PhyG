#!/usr/bin/env python3

import glob, sys,time
import xml.etree.ElementTree as ET

from itertools import product as ip

from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def check_file_exists(file: str) -> bool:
    try:
        filesize = Path(file).stat().st_size

    except FileNotFoundError:
        return False

    if filesize > 0:
        return True

    else:
        return False

def prep_gc_dir(
    out_dir: str,
    taxon_name: str,
    taxon_code: str) -> str:

    q_tax = taxon_name
    if taxon_code:
        q_tax = taxon_code

    if not out_dir:
        gen_code_dir = f'Eval_Genetic_Code_PhyG/'
    else:
        gen_code_dir = f'{out_dir}/Eval_Genetic_Code/'

    Path(gen_code_dir).mkdir(parents = True, exist_ok = True)

    out_xml = f'{gen_code_dir}{q_tax}.DIAMOND_GenCode_Align.xml'

    return gen_code_dir, out_xml


def parse_blast_align(
    start_time,
    out_dir: str,
    out_xml: str,
    fasta_file: str,
    diamond_db: str,
    threads: int = 4,
    verbose: bool = True) -> dict:

    from bin import phyg_orf_related as oc

    xml_prepped = False

    if glob.glob(out_xml) and Path(out_xml).stat().st_size > 0:
        xml_prepped = True

    if not xml_prepped:
        if verbose:
            print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Aligning Query Nucleotides to Proteome Database')

        oc.diamond_og_assign(
            out_dir,
            out_xml,
            fasta_file,
            diamond_db,
            threads,
            eval_gen_code = True)

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Parsing Output XML-Alignment File')

    seq_info = {}

    root = ET.parse(out_xml).getroot()
    dmnd_iters = root.find('BlastOutput_iterations')


    for type_tag in dmnd_iters.findall('Iteration'):
        query_seq = type_tag.find('Iteration_query-def').text

        if len(type_tag[4]) > 0:
            hit_seq = type_tag[4][0][1].text
            qseq_aln = type_tag[4][0][-1][0][-3].text
            hseq_aln = type_tag[4][0][-1][0][-2].text

            qseq_start = int(type_tag[4][0][-1][0][4].text)
            qseq_end = int(type_tag[4][0][-1][0][5].text)
            qframe = int(type_tag[4][0][-1][0][8].text)

            seq_info[query_seq] = [[qseq_start, qseq_end, qframe],[qseq_aln, hseq_aln]]

    return seq_info


def compare_aln(
        qseq: str,
        codon_pairs: list) -> dict:
    n = 0
    a = 0

    aln_cdns = defaultdict(list)

    for i in codon_pairs:
        a += 1
        if '-' not in i:
            aln_cdns[qseq[n*3:n*3+3]].append(i[1])
            n += 1

        elif i[0] != '-' and i[1] == '-':
            n += 1
            # print(a)

    return aln_cdns


def eval_codon_usage(
        fasta_file: str,
        seq_info: dict) -> list:

    total_codons = 0
    final_eval = []

    cdn_eval = {codon:[] for codon in [Seq(''.join(i)) for i in ip(['G','C','A','T'], repeat = 3)]}

    for i in SeqIO.parse(fasta_file,'fasta'):
        if i.id in seq_info:
            orf_info, q_h_aln = seq_info[i.id]

            qseq = i.seq[orf_info[0]-1: orf_info[1]]

            if orf_info[2] < 0:
                qseq = qseq.reverse_complement()

            qseq_cdn_eval = compare_aln(qseq, list(zip(q_h_aln[0], q_h_aln[1])))

            for k, v in qseq_cdn_eval.items():
                if k in cdn_eval:
                    cdn_eval[k] += v
                    total_codons += len(v)

    for k, v in cdn_eval.items():
        check_same = 'Same'

        std_cdn = f'{k.translate()}'
        if std_cdn == '*':
            std_cdn = 'STOP'

        cdn_count = len(v)

        if cdn_count != 0:
            common_aa = max(set(v), key = v.count)
            common_aa_freq = v.count(common_aa)

            if common_aa != std_cdn:
                check_same = 'Reassigned'

            if v.count(common_aa)/cdn_count < 0.33:
                common_aa = common_aa.lower()
            else:
                common_aa = common_aa.upper()

            if cdn_count > total_codons * 0.0001:
                final_eval.append(f'{k}\t{std_cdn}\t{common_aa}\t{check_same}\t{common_aa_freq}\t{cdn_count}\t{total_codons}')
            else:
                if std_cdn == 'STOP':
                    check_same = 'Same'
                    final_eval.append(f'{k}\t{std_cdn}\tSTOP\t{check_same}\t{cdn_count}\t{cdn_count}\t{total_codons}')

        else:
            common_aa = '?'
            common_aa_freq = 0
            if std_cdn == 'STOP':
                final_eval.append(f'{k}\t{std_cdn}\tSTOP\tSame\t{common_aa_freq}\t{cdn_count}\t{total_codons}')
            else:
                final_eval.append(f'{k}\t{std_cdn}\t{common_aa}\tUncertain\t{common_aa_freq}\t{cdn_count}\t{total_codons}')

    final_eval.sort(key=lambda x: (x.split('\t')[1], x.split('\t')[0]))

    return final_eval

    # cdn_count = len(v)
    # aa_freq = sorted(((i, v.count(i)/len(v)) for i in set(v)), key=lambda x: -x[-1])
    # if aa_freq[0][-1] < 0.5:
    #     common_aa = ';'.join([i[0].lower() for i in aa_freq if i[1] > 0.3])
    # else:
    #     common_aa = aa_freq[0][0]
    #
    # if v.count(common_aa)/cdn_count < 0.5:
    #     common_aa = common_aa.lower()
    # else:
    #     common_aa = common_aa.upper()


def eval_genetic_code(
        start_time,
        fasta_file: str,
        diamond_db: str,
        taxon_name: str,
        taxon_code: str = None,
        out_dir: str = None,
        threads: int = 4,
        verbose = True) -> str:

    gen_code_dir, out_xml = prep_gc_dir(
                                out_dir,
                                taxon_name,
                                taxon_code)

    seq_info = parse_blast_align(
                start_time,
                gen_code_dir,
                out_xml,
                fasta_file,
                diamond_db,
                threads,
                verbose)

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Inferring Codons to Amino Acid Translations')

    genetic_code_table = eval_codon_usage(
                            fasta_file,
                            seq_info)
    q_tax = taxon_name
    if taxon_code:
        q_tax = taxon_code

    if verbose:
        print(f'[{timedelta(seconds = round(time.time()-start_time))}]  Saving Inferred Codons to Amino Acid Translations')

    with open(f'{gen_code_dir}{q_tax}.Genetic_Code_Estimate.PhyG.tsv','w+') as w:
        w.write('Codon\tStandard-Amino-Acid\tPredicted-Amino-Acid\tSame-Different\tAA-Frequency\tAA-Count\tTotal-Codons\n')
        w.write('\n'.join(genetic_code_table))

    return out_xml
