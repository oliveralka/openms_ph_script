# openms_ph_script

Tool to extract the location of the phosphorylation site at the peptide and protein level.

Used in combination with OpenMS and the LuciphorAdapter. 

Input:
1) database used for the analysis
2) output of the IdConflictResolver after mzTab export.

Please use as follows:

python extractPhosphoPosition.py --help

Usage: extractPhosphoPosition.py [OPTIONS]

Options:
  -db, --database PATH  Database used for peptide search (.fasta)
  -in, --analysis PATH  TMT10Plex,Phospho data exported from
                        IDConflictResolver and MzTabExporter (.mztab)
  -out, --output PATH   Path to save the output file
  --help                Show this message and exit.

The KNIME Workflow which can be used for TMT, Phospho analysis will be added soon.
