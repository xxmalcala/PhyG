#!/usr/bin/env python3

"""
Controls the flow for PhyG's modules. Currently does NOT resume in-progress
runs. This will be fixed. Additionally, there will be an "adding" data suite of
modules/scripts, that are yet to come [as of 2024-05-11].
"""

import logging,sys, time
from datetime import datetime, timedelta
from pathlib import Path


from bin import phyg_eval_args as eval_args


def check_args():
    args = eval_args.collect_args()

    if (args.transcripts or args.orfs):
        args.add_taxa = True

        if not args.taxon_name:
            print('[ERROR]: Missing a taxon/sample name! Please provide (e.g., "Homo sapiens" or "Me138-q")')
            sys.exit(1)

        if not args.db:
            print('[ERROR]: Missing a gene-family database for gene-family assignment')
            sys.exit(1)

        elif not args.og_delim:
            print('[ERROR]: Missing a delimiter to determine gene-family names in'
                ' from gene-family "hits" report')
            sys.exit(1)

        else:

            args.taxon_name = '_'.join(args.taxon_name).replace(".","").replace('-','_')

            args.add_taxa = True

            args.taxon_dir = f'{args.taxon_name}_PhyG'
            if args.taxon_code:
                args.taxon_dir = f'{args.taxon_code}_PhyG'


            args.remove_isoforms = True

            if args.all_isoforms:
                args.remove_isoforms = False

        # Need a taxon-evluation for "prokaryotic" -- can be mentioned as an option,
        # or checked by the taxonomy script... for now, it's just "off"

            args.prokaryotic = False

            return args

    else:

        args.add_taxa = False

        if args.msa_dir and not args.out_dir:
            args.out_dir = f'{args.msa_dir.rpartition("/")[-1]}_PhyG/'

        if not args.out_dir:
            print('[ERROR]: Missing a project/output directory name! Please provide (e.g., "Phylogeny_Test")')
            sys.exit(1)

        if not args.no_msa_filter:
            args.run_msa = True

            if args.os_kmer:
                args.os_sim = False

            if not args.no_size_filter:
                args.size_filter = True

            if args.size_filter_mean:
                arg.size_filter_median = False

            if args.os_sim_iqr:
                args.os_sim_boots = False

            else:
                if not args.lof:
                    args.os_sim_boots = True

                else:
                    args.os_sim_boots = False

            if not args.msa_dir:
                num_errors = 0
                if not args.prot_dir:
                    print('[ERROR]: Missing a path to peptide directory!')
                    num_errors += 1

                if not args.taxon_list:
                    print('[ERROR]: Missing a list of taxa to include from the peptide directory!')
                    num_errors += 1

                if not args.gf_list:
                    print('[ERROR]: Missing a list of gene families to include from the peptide directory!')
                    num_errors += 1

                if num_errors > 0:
                    sys.exit(1)

        # Maybe make the user pick the approach for MSA refinement -- LOF recommended
        # otherwise, skip all the non-size filtering

        # IF skip-filtering is on, then skip all filter steps

        else:
            args.os_sim_boots = False
            args.size_filter = False

        if args.only_msa:
            args.run_msa = True

        return args


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

    args = check_args()

    # print(args.no_msa_filter)
    # sys.exit()
    # print(args.only_top_hit)
    # sys.exit()

    if args.add_taxa:
        from bin import phyg_add_taxa as add_taxa

        if args.transcripts:
            add_taxa.prep_transcriptomes(
                        args.transcripts,
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
                        args.threads
                        )

        elif args.orfs:
            add_taxa.prep_wgs(
                        args.orfs,
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
                        args.prokaryotic,
                        args.remove_isoforms,
                        args.threads)

    else:
        if args.run_msa:
            from bin import phyg_msa_prep as msa_prep
            # Need a check for the MSA-dir

            aln_dir = msa_prep.prep_gfs_eval_msas(
                        args.out_dir,
                        args.msa_dir,
                        args.prot_dir,
                        args.taxon_list,
                        args.gf_list,
                        args.delim,
                        args.og_delim,
                        args.msa_prog,
                        args.size_filter,
                        args.size_filter_median,
                        args.min_prop,
                        args.max_prop,
                        args.os_guidance,
                        args.os_guidance_bs,
                        args.os_kmer,
                        args.os_kmer_len,
                        args.os_sim,
                        args.os_sim_boots,
                        args.os_sim_iqr,
                        args.os_sim_thresh,
                        args.os_sim_ns,
                        args.os_sim_bs,
                        args.threads)



    # dpnd, args = ea.check_all()
    # print(args)
    # init_log_stats(args)
    #
    # print(f'\n{ea.ascii_logo_vsn()}')

    # prep_msa_files(dpnd['Guidance2'], args)

    print('Thanks for using PhyG!')
