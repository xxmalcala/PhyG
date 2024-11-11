#!/usr/bin/env python3

import sys

from pathlib import Path

from ete3 import NCBITaxa


def search_ncbi_lineage(
        taxon_name: str,
        full_taxonomy: bool = True) -> str:
    """
    Returns NCBI-based taxonomy for a given taxon name

    Parameters
    ----------
    taxon_name:  taxon's name, with white-speces replaced with underscores

    Returns
    ----------
    query_taxonomy:  the NCBI taxonomy for a given taxonomic name
    """
    ncbi = NCBITaxa()

    lineage_final = []
    query_taxonomy = []

    if not Path(f'{Path.home()}/.etetoolkit/taxa.sqlite').is_file():
        ncbi.update_taxonomy_database()

    taxon = ' '.join(taxon_name.split("_")[:2])

    if taxon.split()[1] in ['sp','cf']:
        taxon = taxon.split()[0]

    taxid = ncbi.get_name_translator([taxon])

    if not taxid:
        taxid = ncbi.get_name_translator([taxon.split()[0]])
        if not taxid:
            return 'No NCBI Taxonomy'

    taxon_lineage = ncbi.get_lineage(taxid[taxon][0])

    lineage_names = ncbi.get_taxid_translator(taxon_lineage)

    # Try and keep these major ranks... this is quite variable by taxon...
    ranks_to_check = ['superkingdom','kingdom','clade','phylum','class','order','family','genus']


    for taxid in taxon_lineage:
        if list(ncbi.get_rank([taxid]).values())[0] in ranks_to_check:
            query_taxonomy.append(lineage_names[taxid])

    if full_taxonomy:
        return ';'.join(query_taxonomy)

    # just make a dict of these.... this would be easier to control for now...
    # brief_taxonomy = []
    # if query_taxonomy[2] in ['Metazoa','Fungi']:
    #     domain, supergroup, mjr_clade, mnr_clade = 'Eukaryota', query_taxonomy[1], query_taxonomy[2], query_taxonomy[5]
    #
    # elif query_taxonomy[1] in ['Viridiplantae', 'Rhodophyta', 'Rhodelphea','Glaucocystophyceae']:
    #     domain, supergroup, mjr_clade, mnr_clade = 'Eukaryota', 'Archaeplastida', query_taxonomy[1], query_taxonomy[2]
    #
    # elif 'SAR' in query_taxonomy:
    #     domain, supergroup, mjr_clade, mnr_clade = query_taxonomy[0], query_taxonomy[1], query_taxonomy[2], query_taxonomy[3]
    #
    # else:
    #     domain, supergroup, mjr_clade, mnr_clade = query_taxonomy[0], query_taxonomy[1], query_taxonomy[2]
