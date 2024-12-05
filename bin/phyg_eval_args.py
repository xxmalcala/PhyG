#!/usr/bin/env python3


import argparse, sys

def capture_args():
    # create the top-level parser
    gen_descript = 'Something goes here...\n\nUsage:\n     phyg.py <module> <options>\n\n' \
        'To see module options:\n     phyg.py <module> -h'

    parser = argparse.ArgumentParser(description = gen_descript,
            usage = argparse.SUPPRESS, add_help = False,
            formatter_class = argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--help', '-h', action = "help", help = argparse.SUPPRESS)
    parser.add_argument('--version', action = 'store_true', help = argparse.SUPPRESS)

    # create sub-parsers
    sub_parsers = parser.add_subparsers(title = 'PhyG Modules', description = (
    '''\nadd-taxa        assign gene families to a new dataset\n\n'''
    '''contam          assess putative contamination\n\n'''
    '''msa             multi-sequence alignment methods\n\n'''
    '''tree            phylogenetic tree reconstructions methods\n\n'''
    '''msa-tree        multi-sequence alignment and
                single-gene tree reconstructions\n\n'''
    '''genetic-code    infer genetic code (codon-to-amino acids)\n\n'''),
    help = argparse.SUPPRESS, dest = 'command')

    # create the parser for the "add-taxa" sub-command
    parser_add_taxa = sub_parsers.add_parser('add-taxa', usage = argparse.SUPPRESS, add_help = False,
                        formatter_class = argparse.RawDescriptionHelpFormatter)

    add_taxa = parser_add_taxa.add_argument_group('General Adding-Taxa Options', description = (
    '''--in (-i)             input FASTA file\n\n'''
    '''--out-dir (-o)        output directory name\n\n'''
    '''--taxon-name          taxon name (genus species)\n\n'''
    '''--taxon-code          additional taxonomic code included\n\n'''
    '''--transcripts         FASTA file includes transcripts\n\n'''
    '''--orfs                FASTA file includes CDS/ORFs\n\n'''
    '''--db (-d)             gene-family database (FASTA, DIAMOND, or HMMer format)\n\n'''
    '''--delim               delimiter to use in sequence names
                      (include quotations; default = '|')\n\n'''
    '''--og-delim            delimiter for gene-family names in gene-family database\n\n'''
    '''--gen-code (-g)       translation table to use to translate ORFs (default = 1)\n\n'''
    '''--only-top-hit        keep only the top HSP-scoring hit for each gene-family\n\n'''
    '''--threads (-t)        number of CPU threads to use (default = 4)\n\n'''
    '''--stats               for when you want more stats and figures!\n\n'''
    '''--clean               remove intermediate files\n\n'''
    '''--quiet (-q)          no console output\n\n'''
    '''--gzip (-gz)          tar and gzip TIdeS output\n'''
    '''--help (-h)           show this help message and exit\n\n'''))

    add_orfs = parser_add_taxa.add_argument_group('ORF and CDS Options', description = (
    '''--genbank             use if data are CDSs from GenBank or RefSeq\n'''
    '''--refseq              use if data are CDSs from GenBank or RefSeq\n'''
    '''--all-isoforms        keep all CDS isoforms from a CDS file
                      sourced from GenBank\n'''))

    add_txpts = parser_add_taxa.add_argument_group('Transcriptome Options', description = (
    '''--min-len             minimum transcript length to evaluate\n'''))

    add_blast = parser_add_taxa.add_argument_group('BLAST OG-Assignment Options', description = (
    '''--evalue (-e)         maximum e-value to infer reference ORFs\n\n'''
    '''--subject-cover       minimum alignment length (proportion of the hit)
                      for a "hit" to the gene-family database (default = 60)\n\n'''
    '''--query-cover         minimum alignment length (proportion of the query)
                      for a "hit" to the gene-family database (default = 0)\n\n'''
    '''--min-id              minimum % identity for a "hit" to the gene-family
                      database (default = 0)\n\n'''))

    add_blast = parser_add_taxa.add_argument_group('Hmmer OG-Assignment Options', description = (
        '''--evalue (-e)         maximum e-value to infer reference ORFs\n\n'''))

    add_taxa.add_argument('--help', '-h', action = "help", help = argparse.SUPPRESS)

    add_taxa.add_argument('--input', '--in', '-i', action = 'store', metavar = '[FASTA-file]',
        type = str, help = argparse.SUPPRESS)

    add_taxa.add_argument('--out-dir','-o', action = 'store', metavar = '[output-directory]',
        type = str, help = argparse.SUPPRESS)

    add_taxa.add_argument('--taxon-name', action = 'store', metavar = '[taxon-name]',
        nargs='+', type = str, help = argparse.SUPPRESS)

    add_taxa.add_argument('--taxon-code', action = 'store', metavar = '[taxon-code]',
        type = str, help = argparse.SUPPRESS)

    add_taxa.add_argument('--transcripts', action = 'store_true', help = argparse.SUPPRESS)

    add_taxa.add_argument('--orfs', action = 'store_true', help = argparse.SUPPRESS)

    add_taxa.add_argument('--genbank', action = 'store_true', help = argparse.SUPPRESS)

    add_taxa.add_argument('--refseq', action = 'store_true', help = argparse.SUPPRESS)

    add_taxa.add_argument('--db', '-d', action = 'store', metavar = '[OG Database]',
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

    add_taxa.add_argument('--threads','-t', action = 'store', default = 4,
        metavar = '[Threads]', type = int, help = argparse.SUPPRESS)

    add_taxa.add_argument('--stats', action = 'store_true', help = argparse.SUPPRESS)

    add_taxa.add_argument('--clean', action = 'store_true', help = argparse.SUPPRESS)

    add_taxa.add_argument('--quiet', '-q', action = 'store_false', help = argparse.SUPPRESS)

    add_taxa.add_argument('--gzip','-gz', action = 'store_true', help = argparse.SUPPRESS)


    # create the parser for the "contamination" sub-command
    parser_contam = sub_parsers.add_parser('contam', usage = argparse.SUPPRESS, add_help = False,
                        formatter_class = argparse.RawDescriptionHelpFormatter)

    contam = parser_contam.add_argument_group('General Contamination Options', description = (
    '''--in (-i)             input FASTA file\n\n'''
    '''--out-dir (-o)        output directory name\n\n'''
    '''--taxon-name          taxon name (genus species)\n\n'''
    '''--taxon-code          additional taxonomic code included\n\n'''
    '''--transcripts         FASTA file includes transcripts\n\n'''
    '''--orfs                FASTA file includes predicted CDS/ORFs\n\n'''
    '''--genbank             use if data are CDSs from GenBank or RefSeq\n'''
    '''--refseq              use if data are CDSs from GenBank or RefSeq\n'''
    '''--db (-d)             diagnostic gene-family database\n\n'''
    '''--gen-code (-g)       translation table to use to translate ORFs (default = 1)\n\n'''
    '''--kmer (-k)           k-mer based clustering (default)\n\n'''
    '''--phylo (-p)          phylogenetic based evaluation\n\n'''
    '''--threads (-t)        number of CPU threads to use (default = 4)\n\n'''
    '''--clean               remove intermediate files\n\n'''
    '''--quiet (-q)          no console output\n\n'''
    '''--gzip (-gz)          tar and gzip TIdeS output\n\n'''
    '''--help (-h)           show this help message and exit\n'''))

    contam_kmer = parser_contam.add_argument_group('K-mer Based Options', description = (
    '''--cluster-threshold   distance threshold for calling clusters (default = 3)\n\n'''
    '''--eval-threshold      representative proportion of a clade [e.g. Rhizaria'] to
                      assign a clade-name to a cluster (default = 0.6)\n\n'''))

    contam_phylo = parser_contam.add_argument_group('Phylogenetic Tree Based Options', description = (
    '''--ref-msa (-m)            diagnostic mutli-sequence alignments\n\n'''
    '''--ref-trees (-rt)         reference diagnostic phylogenies\n\n'''))

    contam.add_argument('--help', '-h', action = "help", help = argparse.SUPPRESS)

    contam.add_argument('--input', '--in', '-i', action = 'store', metavar = '[FASTA-file]',
        type = str, help = argparse.SUPPRESS)

    contam.add_argument('--out-dir','-o', action = 'store', metavar = '[output-directory]',
        type = str, help = argparse.SUPPRESS)

    contam.add_argument('--taxon-name', action = 'store', metavar = '[taxon-name]',
        nargs='+', type = str, help = argparse.SUPPRESS)

    contam.add_argument('--taxon-code', action = 'store', metavar = '[taxon-code]',
        type = str, help = argparse.SUPPRESS)

    contam.add_argument('--transcripts', action = 'store_true', help = argparse.SUPPRESS)

    contam.add_argument('--orfs', action = 'store_true', help = argparse.SUPPRESS)

    contam.add_argument('--genbank', action = 'store_true', help = argparse.SUPPRESS)

    contam.add_argument('--refseq', action = 'store_true', help = argparse.SUPPRESS)

    contam.add_argument('--db', '-d', action = 'store', metavar = '[Diagnostic Database]',
        type = str, help = argparse.SUPPRESS)

    contam.add_argument('--gen-code', '-g', action = 'store', default = '1',
        metavar = '[translation-table]', type = str, help = argparse.SUPPRESS)

    contam.add_argument('--kmer', '-k', action = 'store_false', help = argparse.SUPPRESS)

    contam.add_argument('--cluster-threshold', action = 'store', default = 3,
        metavar = '[cluster-thresh]', type = int, help = argparse.SUPPRESS)

    contam.add_argument('--eval-threshold', action = 'store', default = 0.6,
        metavar = '[eval-thresh]', type = int, help = argparse.SUPPRESS)

    contam.add_argument('--phylo', '-p', action = 'store_true', help = argparse.SUPPRESS)

    contam.add_argument('--ref-msa', '-m', action = 'store', metavar = '[diagnostic msas]',
        type = str, help = argparse.SUPPRESS)

    contam.add_argument('--ref-trees', '-rt', action = 'store', metavar = '[diagnostic trees]',
        type = str, help = argparse.SUPPRESS)

    contam.add_argument('--blen-mode', '-bl', action = 'store', default = 'median',
        metavar = '[branch-length-eval]', type = str, help = argparse.SUPPRESS)

    contam.add_argument('--threads','-t', action = 'store', default = 4,
        metavar = '[Threads]', type = int, help = argparse.SUPPRESS)

    contam.add_argument('--clean', action = 'store_true', help = argparse.SUPPRESS)

    contam.add_argument('--quiet', '-q', action = 'store_false', help = argparse.SUPPRESS)

    contam.add_argument('--gzip', '-gz', action = 'store_true', help = argparse.SUPPRESS)


    # create the parse for the "msa" sub-command
    parser_msa = sub_parsers.add_parser('msa', usage = argparse.SUPPRESS, add_help = False,
                        formatter_class = argparse.RawDescriptionHelpFormatter)

    msa_aln = parser_msa.add_argument_group('General Alignment Options', description = (
    '''--out-dir (-o)        output directory name\n\n'''
    '''--msa-dir (-m)        directory of unaligned MSAs\n\n'''
    '''--prot-dir (-p)       directory of FASTA files of assigned gene families
                      per taxon\n\n'''
    '''--gf-list (-g)        list of gene families to align\n\n'''
    '''--taxon-list          list of taxa to include in MSAs\n\n'''
    '''--og-delim (-d)       delimiter for gene-family names\n\n'''
    '''--msa-prog            alignment tool to use (MAFFT or WITCH; default = MAFFT)\n\n'''
    '''--gap-thresh          gaps threshold for trimming column (0 - 1.0; default = 0.9)\n\n'''
    '''--clipkit             use clipkit and "smart-gap" for gap removal (default)\n\n'''
    '''--threads (-t)        number of CPU threads to use (default = 4)\n\n'''
    '''--clean               remove intermediate files\n\n'''
    '''--quiet (-q)          no console output\n\n'''
    '''--gzip (-gz)          tar and gzip TIdeS output\n\n'''
    '''--stats               for when you want more stats and figures!\n\n'''
    '''--help (-h)           show this help message and exit\n'''))

    msa_aln_filt = parser_msa.add_argument_group('General Filtering Options', description = (
    '''--skip-filter         do not remove any sequences from the MSA\n\n'''
    '''--skip-size-filter    do not remove sequences if "too short/long"\n'''))

    msa_size_filt = parser_msa.add_argument_group('Size Filtering Options', description = (
    '''--size-filter-mean    filter sequences based on mean length of the gene-family\n\n'''
    '''--size-filter-median  filter sequences based on median length of the
                      gene-family (default)\n\n'''
    '''--min-prop            minimum proportion of the mean/median gene-family
                      length to keep peptides (default = 0.33)\n\n'''
    '''--max-prop            maximum proportion of the mean/median gene-family
                      length to keep peptides (default = 2.0)\n'''))

    msa_hom_filt = parser_msa.add_argument_group('Homology Assessment and Filtering Options', description = (
    '''--guidance            absolute path to GUIDANCE2 to identify and remove
                      outlier peptide sequences\n\n'''
    '''--guidance-bs         number of bootstraps for GUIDANCE2 outlier removal
                      (default = 10)\n\n'''
    '''--kmer (-k)           k-mer frequencies and unsupervised learning (local
                      outlier factor) to infer outlier peptide sequences\n\n'''
    '''--ngram (-n)          n-gram for outlier sequence removal (default = 5)\n\n'''
    '''--sim (-s)            "sentence similarity" to identify and remove
                      outlier peptide sequences\n\n'''
    '''--sim-thresh          threshold (above which) sequences are removed
                      (from 2 - 10, default = 3)\n\n'''
    '''--sim-iqr            interquartile range of "sentence similarities" as
                      basis for scoring sequences\n\n'''
    '''--sim-bs             number of bootstraps for "sentence similarity" metrics
                      (default = 100)\n\n'''
    '''--sim-ns             number of peptides to sub-sample to create distribution
                      of "sentence similarities" (default = 100)\n\n'''
    '''--lof                use unsupervised learning (local outlier factor) to infer
                     outlier peptides with "sentence similarity"\n\n'''))

    msa_aln.add_argument('--help', '-h', action = "help", help = argparse.SUPPRESS)

    msa_aln.add_argument('--out-dir','-o', action = 'store', metavar = '[output-directory]',
        type = str, help = argparse.SUPPRESS)

    msa_aln.add_argument('--msa-dir', '-m', action = 'store', metavar = '[msa-directory]',
                type = str, help = argparse.SUPPRESS)

    msa_aln.add_argument('--prot-dir', '-p', action = 'store', metavar = '[peptide-directory]',
                type = str, help = argparse.SUPPRESS)

    msa_aln.add_argument('--gf-list', '-g', action = 'store', metavar = '[gene-family-list]',
                type = str, help = argparse.SUPPRESS)

    msa_aln.add_argument('--taxon-list', action = 'store', metavar = '[taxon-list]',
                type = str, help = argparse.SUPPRESS)

    msa_aln.add_argument('--og-delim', '-d', action = 'store', metavar = '[og-delimiter]',
                type = str, help = argparse.SUPPRESS)

    msa_aln.add_argument('--msa-prog', action = 'store', metavar = '[msa-program]',
                type = str, default = 'MAFFT', help = argparse.SUPPRESS)

    msa_aln.add_argument('--gap-thresh', action = 'store', metavar = '[gap-threshold]',
                type = float, default = 0.9, help = argparse.SUPPRESS)

    msa_aln.add_argument('--clipkit', action = 'store_false', help = argparse.SUPPRESS)

    msa_aln.add_argument('--skip-filter', action = 'store_true', help = argparse.SUPPRESS)

    msa_aln.add_argument('--skip-size-filter', action = 'store_true', help = argparse.SUPPRESS)

    msa_aln.add_argument('--size-filter-mean', action = 'store_true', help = argparse.SUPPRESS)

    msa_aln.add_argument('--size-filter-median', action = 'store_false', help = argparse.SUPPRESS)

    msa_aln.add_argument('--min-prop', action = 'store', default = 0.33,
                metavar = '[min-prop-size-filter]', type = float, help = argparse.SUPPRESS)

    msa_aln.add_argument('--max-prop', action = 'store', default = 2.0,
                metavar = '[max-prop-size-filter]', type = float, help = argparse.SUPPRESS)

    msa_aln.add_argument('--guidance', action = 'store_true', help = argparse.SUPPRESS)

    msa_aln.add_argument('--guidance-bs', action = 'store', default = 10,
                metavar = '[guidance2-boot-straps]', type = float, help = argparse.SUPPRESS)

    msa_aln.add_argument('--kmer', '-k', action = 'store_true', help = argparse.SUPPRESS)

    msa_aln.add_argument('--ngram','-n', action = 'store', default = 5,
                metavar = '[filter-kmer-length]', type = float, help = argparse.SUPPRESS)

    msa_aln.add_argument('--sim', action = 'store_true', help = argparse.SUPPRESS)

    msa_aln.add_argument('--sim-thresh', action = 'store', default = 3, help = argparse.SUPPRESS)

    msa_aln.add_argument('--sim-iqr', action = 'store_true', help = argparse.SUPPRESS)

    msa_aln.add_argument('--sim-bs', action = 'store', default = 100,
                metavar = '[filter-sent-sim-boot-strap]', type = float, help = argparse.SUPPRESS)

    msa_aln.add_argument('--sim-ns', action = 'store', default = 100,
                metavar = '[filter-sent-sim-num-seqs]', type = float, help = argparse.SUPPRESS)

    msa_aln.add_argument('--lof', action = 'store_true', help = argparse.SUPPRESS)

    msa_aln.add_argument('--threads','-t', action = 'store', default = 4,
        metavar = '[Threads]', type = int, help = argparse.SUPPRESS)

    msa_aln.add_argument('--clean', action = 'store_true', help = argparse.SUPPRESS)

    msa_aln.add_argument('--quiet', '-q', action = 'store_false', help = argparse.SUPPRESS)

    msa_aln.add_argument('--gzip', '-gz', action = 'store_true', help = argparse.SUPPRESS)

    msa_aln.add_argument('--stats', action = 'store_true', help = argparse.SUPPRESS)


    # create the parse for the "tree" sub-command
    parser_tree = sub_parsers.add_parser('tree', usage = argparse.SUPPRESS, add_help = False,
                        formatter_class = argparse.RawDescriptionHelpFormatter)

    tree_gen = parser_tree.add_argument_group('General Phylogenetic Reconstruction Options', description = (
    '''--msa-dir (-m)        directory of aligned MSAs\n\n'''
    '''--out-dir (-o)        output directory name\n\n'''
    '''--threads (-t)        number of CPU threads to use (default = 4)\n\n'''
    '''--clean               remove intermediate files\n\n'''
    '''--quiet (-q)          no console output\n\n'''
    '''--gzip (-gz)          tar and gzip TIdeS output\n\n'''
    '''--help (-h)           show this help message and exit\n'''))

    tree_phy = parser_tree.add_argument_group('Phylogenetic Reconstruction Options', description = (
    '''--fast-tree           phylogenetic tree construction with FastTree2\n\n'''
    '''--iqtree              phylogenetic tree construction with IQTree2 (default)\n\n'''
    '''--model               specify evolutionary model for phylogenetic tree
                      reconstruction (LG+G, JTT, WAG, etc)\n\n'''
    '''--model-sub           restrict to those AA models designed for either
                      nuclear, mitochondrial, chloroplast or viral
                      (default = nuclear)\n\n'''
    '''--ufb, -B             number of ultrafast bootstraps
                      (default = 1000)\n\n'''
    '''--alrt                number of replicates (>=1000) to perform SH-like
                      approximate likelihood ratio test (SH-aLRT)\n\n'''
    '''--asteroid            species-tree construction with ASTEROID\n\n'''
    '''--aster               species-tree construction with ASTER\n\n'''))


    tree_gen.add_argument('--help', '-h', action = "help", help = argparse.SUPPRESS)

    tree_gen.add_argument('--out-dir','-o', action = 'store', metavar = '[output-directory]',
                type = str, help = argparse.SUPPRESS)

    tree_gen.add_argument('--msa-dir', '-m', action = 'store', metavar = '[msa-directory]',
                type = str, help = argparse.SUPPRESS)

    tree_gen.add_argument('--clean', action = 'store_true', help = argparse.SUPPRESS)

    tree_gen.add_argument('--quiet', '-q', action = 'store_false', help = argparse.SUPPRESS)

    tree_gen.add_argument('--gzip', '-gz', action = 'store_true', help = argparse.SUPPRESS)

    tree_gen.add_argument('--model', action = 'store', default = 'MFP',
        metavar = '[model]', type = str, help = argparse.SUPPRESS)

    tree_gen.add_argument('--model-sub', action = 'store', default = 'nuclear',
        metavar = '[model-set]', type = str, help = argparse.SUPPRESS)

    tree_gen.add_argument('--ufb', '-B', action = 'store', default = 1000,
        metavar = '[ultrafast bootstrap]', type = int, help = argparse.SUPPRESS)

    tree_gen.add_argument('--alrt', action = 'store_true', help = argparse.SUPPRESS)

    tree_gen.add_argument('--asteroid', action = 'store_true', help = argparse.SUPPRESS)

    tree_gen.add_argument('--aster', action = 'store_true', help = argparse.SUPPRESS)

    tree_gen.add_argument('--fast-tree', action = 'store_true', help = argparse.SUPPRESS)

    tree_gen.add_argument('--iqtree', action = 'store_false', help = argparse.SUPPRESS)

    # create the parse for the "msa-tree" sub-command
    parser_msa_tree = sub_parsers.add_parser('msa-tree', usage = argparse.SUPPRESS, add_help = False,
                        formatter_class = argparse.RawDescriptionHelpFormatter)


    msa_tree = parser_msa_tree.add_argument_group('General MSA-Tree Options', description = (
    '''--out-dir (-o)        output directory name\n\n'''
    '''--threads (-t)        number of CPU threads to use (default = 4)\n\n'''
    '''--clean               remove intermediate files\n\n'''
    '''--quiet (-q)          no console output\n\n'''
    '''--gzip (-gz)          tar and gzip TIdeS output\n\n'''
    '''--stats               for when you want all the stats and figures!\n\n'''
    '''--help (-h)           show this help message and exit\n'''))

    msa_tree_aln = parser_msa_tree.add_argument_group('General Alignment Options', description = (
    '''--msa-dir (-m)        directory of unaligned MSAs\n\n'''
    '''--prot-dir (-p)       directory of FASTA files of assigned gene families
                      per taxon\n\n'''
    '''--gf-list (-g)        list of gene families to align\n\n'''
    '''--taxon-list          list of taxa to include in MSAs\n\n'''
    '''--og-delim (-d)       delimiter for gene-family names\n\n'''
    '''--msa-prog            alignment tool to use (MAFFT or MUSCLE; default = MAFFT)\n\n'''
    '''--gap-thresh          gaps threshold for trimming column (0 - 1.0; default = 0.9)\n\n'''
    '''--clipkit             use clipkit and "smart-gap" for gap removal (default)\n\n'''
    '''--skip-filter         do not remove any sequences from the MSA\n\n'''
    '''--skip-size-filter    do not remove sequences if "too short/long"\n'''))

    msa_tree_size_filt = parser_msa_tree.add_argument_group('Size Filtering Options', description = (
    '''--size-filter-mean    filter sequences based on mean length of the gene-family\n\n'''
    '''--size-filter-median  filter sequences based on median length of the
                      gene-family (default)\n\n'''
    '''--min-prop            minimum proportion of the mean/median gene-family
                      length to keep peptides (default = 0.33)\n\n'''
    '''--max-prop            maximum proportion of the mean/median gene-family
                      length to keep peptides (default = 2.0)\n'''))

    msa_tree_hom_filt = parser_msa_tree.add_argument_group('Homology Assessment and Filtering Options', description = (
    '''--guidance            absolute path to GUIDANCE2 to identify and remove
                      outlier peptide sequences\n\n'''
    '''--guidance-bs         number of bootstraps for GUIDANCE2 outlier removal
                      (default = 10)\n\n'''
    '''--kmer (-k)           k-mer frequencies and unsupervised learning (local
                      outlier factor) to infer outlier peptide sequences\n\n'''
    '''--ngram (-n)          n-gram for outlier sequence removal (default = 5)\n\n'''
    '''--sim (-s)            "sentence similarity" to identify and remove
                      outlier peptide sequences\n\n'''
    '''--sim-thresh          threshold (above which) sequences are removed
                      (from 2 - 10, default = 3)\n\n'''
    '''--sim-iqr            interquartile range of "sentence similarities" as
                      basis for scoring sequences\n\n'''
    '''--sim-bs             number of bootstraps for "sentence similarity" metrics
                      (default = 100)\n\n'''
    '''--sim-ns             number of peptides to sub-sample to create distribution
                      of "sentence similarities" (default = 100)\n\n'''
    '''--lof                use unsupervised learning (local outlier factor) to infer
                     outlier peptides with "sentence similarity"\n\n'''))

    msa_tree_phy = parser_msa_tree.add_argument_group('General Phylogenetic Reconstruction Options', description = (
    '''--fast-tree           phylogenetic tree construction with FastTree2\n\n'''
    '''--iqtree              phylogenetic tree construction with IQTree2 (default)\n\n'''
    '''--model               specify evolutionary model for phylogenetic tree reconstruction
                      (LG+G, JTT, WAG, etc)\n\n'''
    '''--model-sub            restrict to those AA models designed for either nuclear,
                      mitochondrial, chloroplast or viral (default = nuclear)\n\n'''
    '''--ufb, -B             number of ultrafast bootstraps
                      (default = 1000)\n\n'''
    '''--alrt                UPDATE\n\n'''
    '''--asteroid            species-tree construction with ASTEROID\n\n'''
    '''--aster               species-tree construction with ASTER\n\n'''))

    msa_tree.add_argument('--help', '-h', action = "help", help = argparse.SUPPRESS)

    msa_tree.add_argument('--out-dir','-o', action = 'store', metavar = '[output-directory]',
                type = str, help = argparse.SUPPRESS)

    msa_tree.add_argument('--msa-dir', '-m', action = 'store', metavar = '[msa-directory]',
                type = str, help = argparse.SUPPRESS)

    msa_tree.add_argument('--prot-dir', '-p', action = 'store', metavar = '[peptide-directory]',
                type = str, help = argparse.SUPPRESS)

    msa_tree.add_argument('--gf-list', '-g', action = 'store', metavar = '[gene-family-list]',
                type = str, help = argparse.SUPPRESS)

    msa_tree.add_argument('--taxon-list', action = 'store', metavar = '[taxon-list]',
                type = str, help = argparse.SUPPRESS)

    msa_tree.add_argument('--og-delim', '-d', action = 'store', metavar = '[og-delimiter]',
                type = str, help = argparse.SUPPRESS)

    msa_tree.add_argument('--msa-prog', action = 'store', metavar = '[msa-program]',
                type = str, default = 'MAFFT', help = argparse.SUPPRESS)

    msa_tree.add_argument('--gap-thresh', action = 'store', metavar = '[gap-threshold]',
                type = float, default = 0.9, help = argparse.SUPPRESS)

    msa_tree.add_argument('--clipkit', action = 'store_false', help = argparse.SUPPRESS)

    msa_tree.add_argument('--skip-filter', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--skip-size-filter', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--size-filter-mean', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--size-filter-median', action = 'store_false', help = argparse.SUPPRESS)

    msa_tree.add_argument('--min-prop', action = 'store', default = 0.33,
                metavar = '[min-prop-size-filter]', type = float, help = argparse.SUPPRESS)

    msa_tree.add_argument('--max-prop', action = 'store', default = 2.0,
                metavar = '[max-prop-size-filter]', type = float, help = argparse.SUPPRESS)

    msa_tree.add_argument('--guidance', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--guidance-bs', action = 'store', default = 10,
                metavar = '[guidance2-boot-straps]', type = float, help = argparse.SUPPRESS)

    msa_tree.add_argument('--kmer', '-k', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--ngram','-n', action = 'store', default = 5,
                metavar = '[filter-kmer-length]', type = float, help = argparse.SUPPRESS)

    msa_tree.add_argument('--sim', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--sim-thresh', action = 'store', default = 3, help = argparse.SUPPRESS)

    msa_tree.add_argument('--sim-iqr', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--sim-bs', action = 'store', default = 100,
                metavar = '[filter-sent-sim-boot-strap]', type = float, help = argparse.SUPPRESS)

    msa_tree.add_argument('--sim-ns', action = 'store', default = 100,
                metavar = '[filter-sent-sim-num-seqs]', type = float, help = argparse.SUPPRESS)

    msa_tree.add_argument('--lof', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--threads','-t', action = 'store', default = 4,
        metavar = '[Threads]', type = int, help = argparse.SUPPRESS)

    msa_tree.add_argument('--clean', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--quiet', '-q', action = 'store_false', help = argparse.SUPPRESS)

    msa_tree.add_argument('--gzip', '-gz', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--stats', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--model', action = 'store', default = 'MFP',
        metavar = '[model]', type = str, help = argparse.SUPPRESS)

    msa_tree.add_argument('--model-sub', action = 'store', default = 'nuclear',
        metavar = '[model-type]', type = str, help = argparse.SUPPRESS)

    msa_tree.add_argument('--ufb', '-B', action = 'store', default = 1000,
        metavar = '[ultrafast bootstrap]', type = int, help = argparse.SUPPRESS)

    msa_tree.add_argument('--alrt', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--asteroid', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--aster', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--fast-tree', action = 'store_true', help = argparse.SUPPRESS)

    msa_tree.add_argument('--iqtree', action = 'store_false', help = argparse.SUPPRESS)


    parser_gen_code = sub_parsers.add_parser('genetic-code', usage = argparse.SUPPRESS, add_help = False,
                        formatter_class = argparse.RawDescriptionHelpFormatter)

    gen_code = parser_gen_code.add_argument_group('General Genetic Code Options', description = (
    '''--in (-i)             input FASTA file\n\n'''
    '''--out-dir (-o)        output directory name\n\n'''
    '''--taxon-name          taxon name (genus species)\n\n'''
    '''--taxon-code          additional taxonomic code included\n\n'''
    '''--db (-d)             gene-family database (FASTA or DIAMOND)\n\n'''
    '''--threads (-t)        number of CPU threads to use (default = 4)\n\n'''
    '''--clean               remove intermediate files\n\n'''
    '''--quiet (-q)          no console output\n\n'''
    '''--gzip (-gz)          tar and gzip TIdeS output\n\n'''
    '''--help (-h)           show this help message and exit\n'''))


    gen_code.add_argument('--help', '-h', action = "help", help = argparse.SUPPRESS)

    gen_code.add_argument('--input', '--in', '-i', action = 'store', metavar = '[FASTA-file]',
        type = str, help = argparse.SUPPRESS)

    gen_code.add_argument('--out-dir','-o', action = 'store', metavar = '[output-directory]',
        type = str, help = argparse.SUPPRESS)

    gen_code.add_argument('--taxon-name', action = 'store', metavar = '[taxon-name]',
        nargs='+', type = str, help = argparse.SUPPRESS)

    gen_code.add_argument('--taxon-code', action = 'store', metavar = '[taxon-code]',
        type = str, help = argparse.SUPPRESS)

    gen_code.add_argument('--db', '-d', action = 'store', metavar = '[OG Database]',
        type = str, help = argparse.SUPPRESS)

    gen_code.add_argument('--threads','-t', action = 'store', default = 4,
        metavar = '[Threads]', type = int, help = argparse.SUPPRESS)

    gen_code.add_argument('--clean', action = 'store_true', help = argparse.SUPPRESS)

    gen_code.add_argument('--quiet', '-q', action = 'store_false', help = argparse.SUPPRESS)

    gen_code.add_argument('--gzip','-gz', action = 'store_true', help = argparse.SUPPRESS)

    if len(sys.argv) == 1:
#     print(ascii_logo_vsn())
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    return args


def double_check_args(args):
    req_args = []

    if args.command in ['add-taxa', 'contam']:
        if not args.input:
            req_args.append('    --in <input-FASTA-file>')

        if not args.taxon_name:
            req_args.append('    --taxon-name <taxon-name>')

        if not args.db:
            req_args.append('    --db <gene-family-database>')

    if args.command == 'contam':
        if args.phylo:
            if not args.ref_msa:
                req_args.append('    --ref-msa <reference-msa-directory>')
            if not args.ref_trees:
                req_args.append('    --ref-msa <reference-tree-directory>')
            if args.blen_mode.lower() not in ['median','average']:
                req_args.append('    --blen-mode <"median" or "average">')

    if args.command in ['msa','msa-tree', 'tree']:
        if not args.out_dir:
            req_args.append('    --out-dir <output-directory>')

    if args.command in ['msa','msa-tree']:
        if not args.msa_dir:
            if not args.prot_dir:
                req_args.append('    --prot-dir <peptide-directory>')

            if not args.gf_list:
                req_args.append('    --gf-list <gene-family-list>')

            if not args.taxon_list:
                req_args.append('    --taxon-list <taxon-list>')

    if args.command == 'tree':
        if not args.msa_dir:
            req_args.append('    --msa-dir <aligned-msa-directory>')

    if args.command == 'genetic-code':
        if not args.input:
            req_args.append('    --in <input-FASTA-file>')

        if not args.db:
            req_args.append('    --db <gene-family-database>')

        if not (args.taxon_name or args.taxon_code):
            req_args.append('    --taxon-name <taxon-name> OR\n' \
            '    --taxon-code <taxon-code>')

    if req_args:
        print('\nERROR: Missing the following required arguments:')
        print("\n".join(req_args) + '\n')
        print(f'\nFor more information and options:\n    phyg.py {args.command} -h\n')
        return False
    else:
        return True
