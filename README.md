# MetaMetaMerge: merge output from metagenomic taxonomic profiling and binning tools

Vitor C. Piro (vitorpiro@gmail.com)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/metametamerge/README.html)

This tool is part of the MetaMeta Pipeline (https://github.com/pirovc/metameta) but can also be used as a standalone tool.

Input:
------

MetaMetaMerge integrates profiling and binning tools. MetaMetaMerge accepts BioBoxes format directly (https://github.com/bioboxes/rfc/tree/master/data-format) or a .tsv file in the following format:

- Profiling: rank, taxon name or taxid, abundance

Example:

    genus   Methanospirillum        0.0029
    genus   Thermus 0.0029
    genus   568394      0.0029
    species Arthrobacter sp. FB24   0.0835
    species 195      0.0582
    species Mycoplasma gallisepticum        0.0536


- Binning: readid, taxon name or taxid, lenght of sequence assigned

Example:

    M2|S1|R140      354     201
    M2|S1|R142      195     201
    M2|S1|R145      457425  201
    M2|S1|R146      562     201
    M2|S1|R147      1245471 201
    M2|S1|R150      354     201


Database profiles:
------------------

Database profiles can be generated based on a list of accession codes and two scripts (acc2tab.bash and tab2count.bash):

    https://github.com/pirovc/metameta/blob/master/scripts/acc2tab.bash
    https://github.com/pirovc/metameta/blob/master/scripts/tab2count.bash

    cat accession_list.txt | xargs --max-procs=12 -I '{}' bash acc2tab.bash '{}' > output_tax
    tab2count.bash output_tax > output_dbprofile

Example:
--------

    ./MetaMetaMerge.py -i binning_out.tsv profile1.tsv profile2.out -d dbprofile1.out dbprofile2.out dbprofile3.out -t 'tool1,tool2,tool3' -c 'b,p,p' -n names.dmp -e nodes.dmp -m merged.dmp -o output_profile.out


Parameters:
-----------
        usage: MetaMetaMerge.py [-h] -i [<input_files> [<input_files> ...]] -d
                                [<database_profiles> [<database_profiles> ...]] -t
                                <tool_identifier> -c <tool_method> -n <names_file> -e
                                <nodes_file> -m <merged_file> [-b <bins>]
                                [-r <cutoff>] [-f <mode>] [-s <ranks>] -o
                                <output_file> [-p <output_type>]
                                [--output-parsed-profiles] [--detailed] [--verbose]
                                [-v]

        MetaMetaMerge by Vitor C. Piro (vitorpiro@gmail.com, http://github.com/pirovc)

        optional arguments:
          -h, --help            show this help message and exit
          -i [<input_files> [<input_files> ...]], --input-files [<input_files> [<input_files> ...]]
                                Input (binning or profiling) files. Bioboxes or tsv
                                format (see README)
          -d [<database_profiles> [<database_profiles> ...]], --database-profiles [<database_profiles> [<database_profiles> ...]]
                                Database profiles on the same order of the input files
                                (see README)
          -t <tool_identifier>, --tool-identifier <tool_identifier>
                                Comma-separated identifiers on the same order of the
                                input files
          -c <tool_method>, --tool-method <tool_method>
                                Comma-separated methods on the same order of the input
                                files (p -> profiling / b -> binning)
          -n <names_file>, --names-file <names_file>
                                names.dmp from the NCBI Taxonomy database
          -e <nodes_file>, --nodes-file <nodes_file>
                                nodes.dmp from the NCBI Taxonomy database
          -m <merged_file>, --merged-file <merged_file>
                                merged.dmp from the NCBI Taxonomy database
          -b <bins>, --bins <bins>
                                Number of bins. Default: 4
          -r <cutoff>, --cutoff <cutoff>
                                Minimum abundance/Maximum results for each taxonomic
                                level (0: off / 0-1: minimum relative abundance / >=1:
                                maximum number of identifications). Default: 0.0001
          -f <mode>, --mode <mode>
                                Result mode (precise, very-precise, linear, sensitive,
                                very-sensitive, no-cutoff). Default: linear
          -s <ranks>, --ranks <ranks>
                                Comma-separated list of ranks to be independently
                                merged (superkingdom,phylum,class,order,family,genus,s
                                pecies,all). Default: species
          -o <output_file>, --output-file <output_file>
                                Output file
          -p <output_type>, --output-type <output_type>
                                Output type (tsv, bioboxes). Default: bioboxes
          --output-parsed-profiles
                                Output parsed and converted profiles for all input
                                files (without cutoff)
          --detailed            Generate an additional detailed output with individual
                                normalized abundances for each tool, where: 0 -> not
                                identified but present in the database, -1 not present
                                in the database.
          --verbose             Verbose output log
          -v, --version         show program's version number and exit
