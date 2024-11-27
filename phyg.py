#!/usr/bin/env python3

"""
Controls the flow for PhyG's modules. Currently does NOT resume in-progress
runs. This will be fixed. Additionally, there will be an "adding" data suite of
modules/scripts, that are yet to come [as of 2024-05-11].
"""

import logging,sys, time
from datetime import datetime, timedelta
from pathlib import Path

from bin import phyg_add_taxa as at
from bin import phyg_contam as ctm
from bin import phyg_eval_args as ea
from bin import phyg_orf_related as oc
from bin import phyg_ncbi_taxonomy as nt

def eval_taxonomy(taxon_name):
    brief_taxonomy = nt.search_ncbi_lineage(
                            taxon_name,
                            True)

    if brief_taxonomy.split(";")[0] in ['Bacteria','Archaea']:
        return True
    else:
        return False


def eval_args(args):
    if not ea.double_check_args(args):
        sys.exit(1)

    else:
        return args

def add_new_taxa(args):
    args.prok = eval_taxonomy(args.taxon_name)
    if args.prok:
        args.gen_code = '11'

    sttime = time.time()
    tax = args.taxon_name

    if args.taxon_code:
        tax = args.taxon_code

    args.taxon_dir = f'{tax}_AddTaxa_PhyG'

    args.blast_based = True

    if args.db.endswith('.hmm'):
        args.blast_based = False

    if not args.transcripts and not args.orfs:
        if (args.refseq or args.genbank):
            args.transcripts = False

        else:
            args.transcripts = at.guess_data_type(
                                    sttime,
                                    args.input,
                                    args.taxon_dir,
                                    args.taxon_name,
                                    args.taxon_code,
                                    args.gen_code,
                                    args.quiet)

        args.orfs = not args.transcripts

    if args.transcripts:
        at.prep_transcriptomes(
                sttime,
                args.input,
                args.taxon_dir,
                args.taxon_name,
                args.taxon_code,
                args.db,
                args.delim,
                args.og_delim,
                args.min_length,
                args.gen_code,
                args.evalue,
                args.min_id,
                args.max_hits,
                args.subject_cover,
                args.query_cover,
                args.only_top_hit,
                args.blast_based,
                args.threads,
                args.quiet)

    else:
        if args.refseq:
            args.genbank = args.refseq

        args.remove_isoforms = not args.all_isoforms

        at.prep_wgs(
            sttime,
            args.input,
            args.taxon_dir,
            args.taxon_name,
            args.taxon_code,
            args.db,
            args.delim,
            args.og_delim,
            args.gen_code,
            args.evalue,
            args.min_id,
            args.max_hits,
            args.subject_cover,
            args.query_cover,
            args.only_top_hit,
            args.blast_based,
            args.genbank,
            args.prok,
            args.remove_isoforms,
            args.threads,
            args.quiet)

def by_catch_assess(args):

    sttime = time.time()

    args.prok = eval_taxonomy(args.taxon_name)
    if args.prok:
        args.gen_code = '11'

    args.kmer = not args.phylo

    tax = args.taxon_name

    if args.taxon_code:
        tax = args.taxon_code

    args.taxon_dir = f'{tax}_Contam_PhyG/'

    if not args.transcripts and not args.orfs:
        if (args.refseq or args.genbank):
            args.transcripts = False

        else:
            args.transcripts = at.guess_data_type(
                                    sttime,
                                    args.input,
                                    args.taxon_dir,
                                    args.taxon_name,
                                    args.taxon_code,
                                    args.gen_code,
                                    args.quiet)

        args.orfs = not args.transcripts

    ctm.eval_contam(
        sttime,
        args.input,
        args.taxon_name,
        args.taxon_dir,
        args.db,
        args.gen_code,
        args.taxon_code,
        args.transcripts,
        args.genbank,
        args.prok,
        args.phylo,
        args.ref_msa,
        args.ref_trees,
        args.blen_mode.lower(),
        args.cluster_threshold,
        args.eval_threshold,
        args.threads,
        args.quiet)

# def log_info(qparams: list, arg_vars: dict) -> None:
#     for k in qparams:
#         x = f'{k}:'.ljust(22)
#         logging.info(f'{x}{str(arg_vars[k.replace("-","_")] or "None")}')
#
# def init_log_stats(args):
#     arg_vars = vars(args)
#
#     curr_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")
#
#     initdir = f'{arg_vars["msa_dir"]}'
#     outdir = f'{arg_vars["project_name"]}_PhyG'
#
#     Path(outdir).mkdir(parents = True, exist_ok = True)
#
#     outlog = f'{outdir}/tides.log'
#
#     logging.basicConfig(filename = outlog,
#         filemode = 'w+',
#         format = '%(message)s',
#         level = logging.DEBUG)
#
#     logging.info(
#         f'+---------------------------------+\n'\
#         f'|         PhyG Parameters         |\n'\
#         f'+---------------------------------+'
#         )
#
#     logging.info(f'\nRun-Start:            {curr_time}')
#     logging.info(f'Query-Dir:            {initdir.rstrip("/")}')
#     logging.info(f'Output-Dir:           {outdir}')
#
#     if arg_vars['no_filtering']:
#         qparams = ['threads']
#         log_info(qparams, arg_vars)
#         logging.info(f'Homology-Filtering:       NONE')
#
#     elif arg_vars['guidance']:
#         qparams = ['max-iter','seq-score','threads']
#         log_info(qparams, arg_vars)
#         logging.info(f'Homology-Filtering:   Guidance2')
#
#     elif arg_vars['iqr']:
#         qparams = ['threads']
#         log_info(qparams, arg_vars)
#         logging.info(f'Homology-Filtering:   P-dist inter-quartile range')
#
#     else:
#         qparams = ['threads']
#         log_info(qparams, arg_vars)
#         logging.info(f'Homology-Filtering:   Bootstrap p-dist')
#
#     return outlog


if __name__ == '__main__':


    # print('Currently re-working argument handling...')
    args = eval_args(ea.capture_args())

    if args.taxon_name:
        args.taxon_name = '_'.join(args.taxon_name).replace(".","")

    if args.command == 'add-taxa':
        if args.quiet:
            print('\n#-------------------------------------------#')
            print('#------------ Adding New Taxon -------------#')
            print('#-------------------------------------------#')

        add_new_taxa(args)

    if args.command == 'contam':
        if args.quiet:
            print('\n#-------------------------------------------#')
            print('#----------- Assessing By-Catch ------------#')
            print('#-------------------------------------------#')

        by_catch_assess(args)

    # print(args)

    if args.quiet:
        print('\nThanks for using PhyG!')
