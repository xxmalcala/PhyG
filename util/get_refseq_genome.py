#!/usr/bin/env python3

#__author__: Xyrus X. Maurer-Alcala
#__last_updated__: 2022-10-20

"""This script will download all CDSs from the representative genome for a user-given taxon.
Additionally, outputs a table with the summary information provided.

Additionally, this script can handle a list of taxon names from a text file.
Under this instance, a SINGLE table will be output along with all the CDSs from the
reference genome assembly."""

import os, sys

from Bio import Entrez as ez
from collections import defaultdict

def get_taxid_num(taxon_name: str) -> int:
    return ez.read(ez.esearch(db = "assembly", term = f'"{taxon_name}"[Organism]'))['IdList']


def get_taxonomy(taxid):
    tmp = ez.efetch(db = 'taxonomy', id = taxid, retmode = 'xml')
    ez_txnmy = ez.read(tmp, validate = True)[0]['Lineage']
    tmp.close()

    return ez_txnmy


def get_genome_info(taxon_name: str) -> dict:
    info_dict = defaultdict(list)

    for tax_num in get_taxid_num(taxon_name):
        tmp_acc_info = ez.esummary(ret_max = None, db = 'assembly', id = tax_num, retmode = 'text')
        ez_info = ez.read(tmp_acc_info, validate = True)
        tmp_acc_info.close()

        for i in ez_info['DocumentSummarySet']['DocumentSummary']:
            if i['RefSeq_category'] == 'representative genome':
                org_name = i['Organism']

                taxid = i['Taxid']
                txnmy = get_taxonomy(taxid)

                acc = i['Synonym']['RefSeq']
                source = 'RefSeq'
                acc_ftp = i['FtpPath_RefSeq']

                if not acc:
                    acc = i['Synonym']['Genbank']
                    if acc:
                        source = 'GenBank'
                        acc_ftp = i['FtpPath_GenBank']
                    else:
                        continue

                info_dict[' '.join(org_name.split()[:2])].append([org_name, taxid, source, acc, f'{acc_ftp}/', txnmy])

    return info_dict


def parse_taxon_file(taxon_file):
    tmp_dict = {}
    for line in open(taxon_file).readlines():
        tmp_dict.update(get_genome_info(line.rstrip()))

    return tmp_dict


def store_table(rs_dict: dict) -> None:
    if os.path.isfile('Genome_Accession_Info.tsv'):
        print('Warning: "Genome_Accession_Info.tsv" already exists...')

        ow = input('Overwrite [y/n]?   ')

        if ow.lower() not in ['y', 'yes']:
            print('\nPlease backup the "Genome_Accession_Info.tsv" file, then run '
                'this script again.\nThis will be less of an issue in the future...\n')

            sys.exit(1)

    with open('Genome_Accession_Info.tsv','w+') as w:
        w.write('Taxon\tFull_Name\tTaxID\tSource\tAccession\tftp_Accession\tNCBI_Taxonomy\n')
        for k, v in rs_dict.items():
            for i in v:
                w.write(f'{k}\t' + '\t'.join(i)+'\n')


if __name__ == '__main__':
    # replace sys with argparse for nargs?
    if len(sys.argv[1:]) == 2:
        taxon_name = sys.argv[1]
        ez.email = sys.argv[2]

    else:
        print('\nUsage:\n    python3 get_refseq_genome.py taxon-name email\n')
        sys.exit(1)

    if taxon_name.endswith('txt'):
        store_table(parse_taxon_file(taxon_name))

    else:
        store_table(get_genome_info(taxon_name))

    os.remove('esummary_assembly.dtd')
