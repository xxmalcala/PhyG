#!/usr/bin/env python3

from pathlib import Path

from ete3 import NCBITaxa


def search_ncbi_lineage(
        taxon_name: str,
        full_taxonomy: bool = True) -> str:
    """
    Returns NCBI-based taxonomy for a given taxon name

    Parameters
    ----------
    taxon_name:    taxon's name, with white-speces replaced with underscores
    full_taxonomy: return taxonomy with all major taxonomic ranks

    Returns
    ----------
    query_taxonomy:    NCBI taxonomy for a given taxonomic name
    reduced_taxonomy:  reduced NCBI taxonomy
    """

    ncbi = NCBITaxa()
    lineage_final = []
    query_taxonomy = []

    if not Path(f'{Path.home()}/.etetoolkit/taxa.sqlite').is_file():
        ncbi.update_taxonomy_database()

    taxon = ' '.join(taxon_name.split("_")[:2])
    taxid = ncbi.get_name_translator([taxon])

    if not taxid:
        if taxon.split()[1].rstrip('.') in ['sp','cf']:
            taxon = taxon.split()[0]
            taxid = ncbi.get_name_translator([taxon])

        else:
            taxon = taxon.split()[0]
            taxid = ncbi.get_name_translator([taxon])

        if not taxid:
            return 'Unknown'

    taxon_lineage = ncbi.get_lineage(taxid[taxon][0])
    lineage_names = ncbi.get_taxid_translator(taxon_lineage)

    # Try and keep these major ranks... this is quite variable by taxon...
    ranks_to_check = ['superkingdom','kingdom','clade','phylum','class','order','family','genus']
    for taxid in taxon_lineage:
        if list(ncbi.get_rank([taxid]).values())[0] in ranks_to_check:
            query_taxonomy.append(lineage_names[taxid])

    if full_taxonomy:
        return ';'.join(query_taxonomy)

    else:
        return reduce_taxonomic_ranks(query_taxonomy)


def reduce_taxonomic_ranks(query_taxonomy: str) -> list:
    """
    Reduces NCBI taxonomy ranks for tree-walking/contamination

    Parameters
    ----------
    query_taxonomy:  list of taxonomic rank names

    Returns
    ----------
    domain:      taxonomic domain
    supergroup:  "Supergroup" or deep taxonomic rank
    mjr_clade:   "Major" taxonomic clade
    mnr_clade:   "Minor" taxonomic clade
    """

    # While some "supergroups" are well supported in some analyses, some are contentious.
    # Contentious supergroups, or those taxa without supergroups, are labeled as "Orphan"
    # for simplicity.
    useful_major_ranks = {'Ancyromonadida':['Orphan','Ancyromonadida'], 'Apusozoa':['Orphan','Apusozoa'],
        'Breviata':['Orphan','Breviata'], 'CRuMs':['Orphan','CRuMs'], 'Cryptophyceae':['Cryptista','Cryptophyceae'],
        'Meteora':['Orphan','Orphan incertae sedis'], 'Telonemida':['Orphan','Telonemida'],
        'Palpitomonas':['Cryptista','Cryptista incertae sedis'], 'Picozoa':['Orphan','Orphan incertae sedis'],
        'Microheliella':['Orphan','Orphan incertae sedis']}

    # Handles eukaryotes
    if query_taxonomy[0] == 'Eukaryota':
        if 'Opisthokonta' in query_taxonomy:
            if query_taxonomy[2] in ['Metazoa']:
                try:
                    domain, supergroup, mjr_clade, mnr_clade = query_taxonomy[0], query_taxonomy[1], query_taxonomy[2], query_taxonomy[5]
                except IndexError:
                    domain, supergroup, mjr_clade, mnr_clade = query_taxonomy[0], query_taxonomy[1], query_taxonomy[2], query_taxonomy[-1]
            else:
                domain, supergroup, mjr_clade, mnr_clade = query_taxonomy[0], query_taxonomy[1], query_taxonomy[2], query_taxonomy[3]

        elif query_taxonomy[1] in ['Viridiplantae', 'Rhodophyta', 'Rhodelphea','Glaucocystophyceae']:
            if len(query_taxonomy) > 2:
                domain, supergroup, mjr_clade, mnr_clade = query_taxonomy[0], 'Archaeplastida', query_taxonomy[1], query_taxonomy[2]

            else:
                domain, supergroup, mjr_clade, mnr_clade = query_taxonomy[0], 'Archaeplastida', query_taxonomy[1], f'{query_taxonomy[1]} incertae sedis'

        elif 'Sar' in query_taxonomy:
            domain, supergroup, mjr_clade, mnr_clade = query_taxonomy[0], query_taxonomy[1].upper(), query_taxonomy[2].replace("Stramenopiles","Stramenopila"), query_taxonomy[3]

        elif query_taxonomy[1] in useful_major_ranks:
            domain = 'Eukaryota'
            mnr_clade = query_taxonomy[-1]
            supergroup, mjr_clade = useful_major_ranks[query_taxonomy[1]]

        elif len(query_taxonomy) > 3:
            domain, supergroup, mjr_clade, mnr_clade = query_taxonomy[0], query_taxonomy[1], query_taxonomy[2], query_taxonomy[3]

        else:
            domain, supergroup, mjr_clade, mnr_clade = 'Unknown', 'Unknown', 'Unknown', 'Unknown'

    # Simple handling for prokarytic datasets
    elif query_taxonomy[0] in ['Bacteria', 'Archaea']:
        if len(query_taxonomy) > 3:
            domain, supergroup, mjr_clade, mnr_clade = query_taxonomy[0], query_taxonomy[1], query_taxonomy[2], query_taxonomy[3]

        elif len(query_taxonomy) == 3:
            domain, supergroup, mjr_clade, mnr_clade = query_taxonomy[0], query_taxonomy[1], query_taxonomy[2], query_taxonomy[2]

        elif len(query_taxonomy) == 2:
            domain, supergroup, mjr_clade, mnr_clade = query_taxonomy[0], f'{query_taxonomy[0]} incertae sedis', f'{query_taxonomy[0]} Incertae sedis', query_taxonomy[1]

        else:
            domain, supergroup, mjr_clade, mnr_clade = 'Unknown', 'Unknown', 'Unknown', 'Unknown'

    else:
        domain, supergroup, mjr_clade, mnr_clade = 'Unknown', 'Unknown', 'Unknown', 'Unknown'

    return domain, supergroup, mjr_clade, mnr_clade
