#!/usr/bin/env python3

import sys

from ete3 import NCBITaxa


def search_ncbi_lineage(taxon_name: str) -> list:
    """
    Returns NCBI-based taxonomy for a given taxon name

    Parameters
    ----------
    taxon_name:  taxon's name, with white-speces replaced with underscores

    Returns
    ----------
    query_taxonomy:  taxonomy as a list, with each rank
    """

    ncbi = NCBITaxa()

    taxon = ' '.join(taxon_name.split("_")[:2])

    taxid = ncbi.get_name_translator([taxon])

    taxon_lineage = ncbi.get_lineage(taxid[taxon][0])

    # Try and keep these major ranks... this is quite variable by taxon...
    ranks_to_keep = ['superkingdom','clade','phylum','class','order','family','genus']

    lineage_ranks = {v:k for k, v in ncbi.get_rank(taxon_lineage).items() if v in ranks_to_keep}

    lineage_names = ncbi.get_taxid_translator(taxon_lineage)

    query_taxonomy = [lineage_names[lineage_ranks[i]] for i in lineage_ranks]

    return query_taxonomy
