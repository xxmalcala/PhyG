#!/usr/bin/env python3


"""
Basic argument handling for PhyG.

Also ensures that all dependencies are available prior to any given run.
"""

import argparse, logging, shutil, sys
from pathlib import Path


def check_dependencies():
    dependencies = {}
    try:
        from Bio import SeqIO
        dependencies['BioPython'] = 'Check'
    except:
        dependencies['BioPython'] = None
    dependencies['Guidance2'] = shutil.which('guidance.pl')
    dependencies['trimAl'] = shutil.which('trimal')
    dependencies['MAFFT'] = shutil.which('mafft')
    dependencies['CD-HIT'] = shutil.which('cd-hit')
    return dependencies


def collect_args():

    parser = argparse.ArgumentParser(description = f'{ascii_logo_vsn()}\n{usage_msg()}',
            usage=argparse.SUPPRESS, add_help = False,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    g = parser.add_argument_group('Phylogenomic Options', description = (
    '''--msa-dir (-m)        directory of MSAs to run with PhyG\n'''
    '''--project-name (-n)   prefix for the output directory\n'''
    '''--threads (-p)        number of CPU threads (default = 1)\n'''
    '''--guid-param (-gp)    which approach to use for Guidance2 (strict, dynamic, default, quick)\n'''
    '''--max-iter (-max)     maximum number of Guidance2 iterations ("0" would run until no sequences are removed, default = 10)\n'''
    '''--seqscore (-ssc)     minimum sequence score cutoff (between 0.0 - 1.0)\n'''
    # '''--no-guid             runs mafft without MSA processing by guidance\n'''
    # '''--resume (-r)         resume a PhyG run from the last checkpoint\n'''
    # '''--keep-all (-k)       keep all intermediate outputs (note this can be large!)\n'''
    # '''--quiet (-q)          no console output\n'''
    # '''--gzip (-gz)          tar and gzip PhyG output\n'''
    '''--help (-h)           show this help message and exit'''))

    g.add_argument('--help', '-h', action="help", help=argparse.SUPPRESS)

    g.add_argument('--msa-dir', '-m', action = 'store',
        metavar = ('[MSA-Directory]'), type = str, required = True,
        help = argparse.SUPPRESS)

    g.add_argument('--project-name', '-n', action = 'store',
        metavar = ('[Project-Name]'), type = str, required = True,
        help = argparse.SUPPRESS)

    g.add_argument('--threads','-p', action = 'store',
        default = 1, metavar = '[Threads]', type = int,
        help = argparse.SUPPRESS)

    g.add_argument('--max-iter', action = 'store',
        metavar = '[Guidance2-Iterations]', type = int, default = 10,
        help = argparse.SUPPRESS)

    g.add_argument('--seqscore','-ssc', action = 'store',
        metavar = '[Guidance2-Params]', type = str, default = 0.6,
        help = argparse.SUPPRESS)

    g.add_argument('--colscore','-csc', action = 'store',
        metavar = '[Guidance2-Params]', type = float, default = 0.7,
        help = argparse.SUPPRESS)

    g.add_argument('--no-guid', action = 'store_true',
        help = argparse.SUPPRESS)

    g.add_argument('--resume','-r', action = 'store_true',
        help = argparse.SUPPRESS)

    g.add_argument('--keep-all','-k', action = 'store_true',
        help = argparse.SUPPRESS)

    g.add_argument('--quiet','-q', action = 'store_true',
        help = argparse.SUPPRESS)

    g.add_argument('--gzip','-gz', action = 'store_true',
        help = argparse.SUPPRESS)


    # Ensures that just script description provided if no arguments provided
    if len(sys.argv[1:]) == 0:
        print(ascii_logo_vsn())
        print(detail_msg())
        print(f'{usage_msg()}\n')
        sys.exit()

    args = parser.parse_args()

    return args


def ascii_logo_vsn():
    alv_msg = """      ___ _         ___
     | _ \ |_ _  _ / __|
     |  _/ ' \ || | (_ |
     |_| |_||_\_, |\___|
              |__/

     Version 0.1.0
    """
    return alv_msg


def detail_msg():
    dtl_msg = """     PhyG ("Fig") is a modular phylogenomic pipeline to
      ease the generation of "filtered" multi-sequence
      alignments and subsequent phylogenies.\n"""
    return dtl_msg


def usage_msg():
    return """Usage:\n    phyg.py [options] --msa-dir [MSA-Directory]"""

def check_all():
    dpnd = check_dependencies()
    args = collect_args()
    # print(dpnd)
    # print(args)
    return dpnd, args

if __name__ == '__main__':
    dpnd = check_dependencies()
    args = collect_args()
    print(dpnd)
    print(args)
