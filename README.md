# openms_ph_script

Script to extract the location of the Phosphorylation at the Peptide and Protein level.

Input:
1) database used for analysis
2) output of ConsensusMapNormalizer

Please use as follows:

python extractPhosphoPosition.py --help

Usage: extractPhosphoPosition.py [OPTIONS]

Options:
  -db, --database PATH  Database used for Peptide search (.fasta)
  -in, --analysis PATH  TMT10Plex,Phospho data exported from
                        ConsensusMapNormalizer (.csv)
  -out, --output PATH   Path to save the output file
  --help                Show this message and exit.

The KNIME Workflow which can be used for TMT, Phospho analysis will be added soon.
