#!/usr/bin/env python3


import argparse, sys

def collect_args():

    parser = argparse.ArgumentParser(description = 'Something goes here',
            usage=argparse.SUPPRESS, add_help = False,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    g = parser.add_argument_group('General Options', description = (
    '''--out-dir (-o)        output directory name\n'''
    '''--threads (-t)        number of CPU threads to use (default = 4)\n'''
    '''--clean               remove intermediate files\n'''
    '''--quiet (-q)          no console output\n'''
    '''--gzip (-gz)          tar and gzip TIdeS output\n'''
    '''--help (-h)           show this help message and exit'''))

    g.add_argument('--help', '-h', action="help", help = argparse.SUPPRESS)

    g.add_argument('--out-dir','-o', action = 'store', metavar = '[output-directory]',
        type = str, help = argparse.SUPPRESS)

    g.add_argument('--threads','-t', action = 'store', default = 4,
        metavar = '[Threads]', type = int, help = argparse.SUPPRESS)

    g.add_argument('--quiet','-q', action = 'store_true',
        help = argparse.SUPPRESS)

    g.add_argument('--gzip','-gz', action = 'store_true',
        help = argparse.SUPPRESS)

    g.add_argument('--version', action = 'store_true',
        help = argparse.SUPPRESS)

    add_taxa = parser.add_argument_group('Adding-Taxa Options', description = (
    '''--taxon-name          taxon name (genus species) [REQUIRED]\n'''
    '''--taxon-code          additional taxonomic code include\n'''
    '''--transcripts         assign gene families to FASTA file of transcripts\n'''
    '''--orfs                assign gene families to FASTA file of CDS/ORFs\n'''
    '''--genbank             use if data are CDSs from GenBank or RefSeq\n'''
    '''--refseq              use if data are CDSs from GenBank or RefSeq\n'''
    '''--all-isoforms        keep all CDS isoforms from a CDS file sourced from GenBank\n'''
    '''--db (-db)            gene-family database (FASTA, DIAMOND, or HMMer format)\n'''
    '''--min-len             minimum transcript length to evaluate\n'''
    '''--evalue (-e)         maximum e-value to infer reference ORFs\n'''
    '''--subject-cover       minimum alignment length (proportion of the hit)
                      for a "hit" to the gene-family database (default = 60)\n'''
    '''--query-cover         minimum alignment length (proportion of the query)
                      for a "hit" to the gene-family database (default = 0)\n'''
    '''--min-id              minimum % identity for a "hit" to the gene-family
                      database (default = 0)\n'''
    '''--delim               delimiter to use in sequence names
                      (include quotations; default = '|')\n'''
    '''--og-delim            delimiter for gene-family names in gene-family database\n'''
    '''--gen-code (-g)       translation table to use to translate ORFs (default = 1)\n'''
    '''--only-top-hit        keep only the top HSP-scoring hit for each gene-family\n\n'''))


    add_taxa.add_argument('--taxon-name', action = 'store', metavar = '[taxon-name]',
                nargs='+', type = str, help = argparse.SUPPRESS)

    add_taxa.add_argument('--taxon-code', action = 'store', metavar = '[taxon-code]',
                type = str, help = argparse.SUPPRESS)

    add_taxa.add_argument('--transcripts', action = 'store', help = argparse.SUPPRESS)

    add_taxa.add_argument('--orfs', action = 'store', help = argparse.SUPPRESS)

    add_taxa.add_argument('--genbank', action = 'store_true', help = argparse.SUPPRESS)

    add_taxa.add_argument('--refseq', action = 'store_true', help = argparse.SUPPRESS)

    add_taxa.add_argument('--db', '-db', action = 'store', metavar = '[OG Database]',
                type = str, help = argparse.SUPPRESS)

    add_taxa.add_argument('--evalue', '-e', action = 'store', default = 1e-10,
                metavar = '[e-value]', type = float, help = argparse.SUPPRESS)

    add_taxa.add_argument('--subject-cover', action = 'store', default = 60,
                metavar = '[subject-cover]', type = float, help = argparse.SUPPRESS)

    add_taxa.add_argument('--query-cover', action = 'store', default = 0,
                metavar = '[query-cover]', type = float, help = argparse.SUPPRESS)

    add_taxa.add_argument('--min-id', action = 'store', default = 0,
                metavar = '[perc-identity]', type = float, help = argparse.SUPPRESS)

    add_taxa.add_argument('--min-length', action = 'store', default = 200,
                metavar = '[min-transcript-length]', type = int, help = argparse.SUPPRESS)

    add_taxa.add_argument('--delim', action = 'store', default = '|',
                metavar = '[name-delimiter]', type = str, help = argparse.SUPPRESS)

    add_taxa.add_argument('--og-delim', action = 'store', default = '|',
                metavar = '[gene-family-name-delimiter]', type = str, help = argparse.SUPPRESS)

    add_taxa.add_argument('--gen-code', '-g', action = 'store', default = '1',
                metavar = '[translation-table]', type = str, help = argparse.SUPPRESS)

    add_taxa.add_argument('--max-hits', action = 'store', default = 1,
                metavar = '[max-gf-hits]', type = int, help = argparse.SUPPRESS)

    add_taxa.add_argument('--only-top-hit', action = 'store_true', help = argparse.SUPPRESS)

    add_taxa.add_argument('--all-isoforms', action = 'store_true', help = argparse.SUPPRESS)

    mng_aln = parser.add_argument_group('Alignment Options', description = (
    '''--only-msa            runs just the multi-sequence alignment steps\n'''
    '''--msa-dir             directory of unaligned MSAs\n'''
    '''--prot-dir            directory of FASTA files of assigned gene families
                      per taxon [REQUIRED]\n'''
    '''--gf-list             list of gene families to align [REQUIRED]\n'''
    '''--taxon-list          list of taxa to include in MSAs\n'''
    '''--og-delim            delimiter for gene-family names for grabbing peptides\n'''
    '''--no-size-filter      do not remove sequences if "too short/long"\n'''
    '''--no-msa-filter       do not remove any sequences from the MSA\n'''
    '''--size-filter-mean    filter sequences based on mean length of the gene-family\n'''
    '''--size-filter-median  filter sequences based on median lengths of
                      the gene-family (default)\n'''
    '''--min-prop            minimum proportion of the mean/median gene-family
                      length to keep peptides (default = 0.33)\n'''
    '''--max-prop            maximum proportion of the mean/median gene-family
                      length to keep peptides (default = 2.0)\n'''
    '''--os-guidance         use GUIDANCE2 to identify and remove outlier peptides
                      from the MSA\n'''
    '''--os-guidance-bs      number of bootstraps for GUIDANCE2 outlier removal
                      (default = 10)\n'''
    '''--os-kmer             use k-mer based NLP+ML approach to identify and
                      remove outlier sequences\n'''
    '''--os-kmer-len         size of k-mer to use for k-mer based outlier sequence
                      removal (default = 5)\n'''
    '''--os-sim              use "sentence similarity" NLP approach to identify
                      and remove outlier sequences (default)\n'''
    '''--os-sim-thresh       threshold (above which) sequences are removed
                      (from 2 - 10, default = 3)\n'''
    '''--os-sim-iqr          use interquartile range of sentence similarities as
                      basis for scoring sequences\n'''
    '''--os-sim-bs           number of bootstraps for sentence-similarity approach
                      (default = 100)\n'''
    '''--os-sim-ns           number of peptides to sub-sample to create distribution
                      of sentence similarities (default = 100)\n'''
    '''--lof                 use unsupervised learning (local outlier factor) to infer
                      outlier peptides (default for k-mer based approaches)\n'''
    '''--msa-prog            alignment tool to use (MAFFT or WITCH; default = MAFFT)\n'''
    ))


    mng_aln.add_argument('--only-msa', action = 'store_true', help = argparse.SUPPRESS)

    mng_aln.add_argument('--msa-dir', action = 'store', metavar = '[msa-directory]',
                type = str, help = argparse.SUPPRESS)

    mng_aln.add_argument('--prot-dir', action = 'store', metavar = '[peptide-directory]',
                type = str, help = argparse.SUPPRESS)

    mng_aln.add_argument('--gf-list', action = 'store', metavar = '[gene-family-list]',
                type = str, help = argparse.SUPPRESS)

    mng_aln.add_argument('--taxon-list', action = 'store', metavar = '[taxon-list]',
                type = str, help = argparse.SUPPRESS)

    mng_aln.add_argument('--no-size-filter', action = 'store_true', help = argparse.SUPPRESS)

    mng_aln.add_argument('--no-msa-filter', action = 'store_true', help = argparse.SUPPRESS)

    mng_aln.add_argument('--size-filter-mean', action = 'store_true', help = argparse.SUPPRESS)

    mng_aln.add_argument('--size-filter-median', action = 'store_false', help = argparse.SUPPRESS)

    mng_aln.add_argument('--min-prop', action = 'store', default = 0.33,
                metavar = '[min-prop-size-filter]', type = float, help = argparse.SUPPRESS)

    mng_aln.add_argument('--max-prop', action = 'store', default = 2.0,
                metavar = '[max-prop-size-filter]', type = float, help = argparse.SUPPRESS)

    mng_aln.add_argument('--os-guidance', action = 'store_true', help = argparse.SUPPRESS)

    mng_aln.add_argument('--os-guidance-bs', action = 'store', default = 10,
                metavar = '[guidance2-boot-straps]', type = float, help = argparse.SUPPRESS)

    mng_aln.add_argument('--os-kmer', action = 'store_true', help = argparse.SUPPRESS)

    mng_aln.add_argument('--os-kmer-len', action = 'store', default = 5,
                metavar = '[filter-kmer-length]', type = float, help = argparse.SUPPRESS)

    mng_aln.add_argument('--os-sim', action = 'store_true', help = argparse.SUPPRESS)

    mng_aln.add_argument('--os-sim_thresh', action = 'store', default = 3, help = argparse.SUPPRESS)

    mng_aln.add_argument('--os-sim-iqr', action = 'store_true', help = argparse.SUPPRESS)

    mng_aln.add_argument('--os-sim-bs', action = 'store', default = 100,
                metavar = '[filter-sent-sim-boot-strap]', type = float, help = argparse.SUPPRESS)

    mng_aln.add_argument('--os-sim-ns', action = 'store', default = 100,
                metavar = '[filter-sent-sim-num-seqs]', type = float, help = argparse.SUPPRESS)

    mng_aln.add_argument('--lof', action = 'store_true', help = argparse.SUPPRESS)

    mng_aln.add_argument('--msa-prog', action = 'store', metavar = '[msa-program]',
                type = str, default = 'MAFFT', help = argparse.SUPPRESS)

    args = parser.parse_args()

    if args.refseq == True:
        args.genbank = True

    args.blast_based = True

    if args.db and args.db.endswith('.hmm'):
        args.blast_based = False


    if len(sys.argv[1:]) == 0:
        print('Custom message coming!')
        # print(ascii_logo_vsn())
        # print(detail_msg())
        # print(f'{usage_msg()}\n')
        sys.exit()

    return args


if __name__ == '__main__':
    print(collect_args())
