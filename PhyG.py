#!/usr/bin/env python3

#__author__: Xyrus X. Maurer-Alcala
#__last_updated__: 2023-04-10

"""Controls the flow for PhyG's modules. Currently does NOT resume in-progress
runs. This will be fixed. Additionally, there will be an "adding" data suite of
modules/scripts, that are yet to come [as of 2023-04-11]."""

import sys
from bin import eval_args as ea
from bin import filt_msa as fm


def prep_msa_files(guid_path, args: list) -> None:

    fin_dir = fm.make_msa_files(guid_path,
                        args.msa_dir,
                        args.project_name,
                        args.max_iter,
                        args.threads,
                        args.seqscore,
                        args.colscore,
                        args.keep_all,
                        args.quiet)

if __name__ == '__main__':
    dpnd, args = ea.check_all()
    print(f'\n{ea.ascii_logo_vsn()}')
    prep_msa_files(dpnd['Guidance2'], args)
