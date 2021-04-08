# miR-Search

miR-Search is a programme developed to identify all 6mer, 7mer and 8mer targets of a given micro-RNA as well as provoiding an interface to download and clean 3'UTR data from the ensembl database to be used in the programme

This project consists of two modules pymart and miR-Search:

pymart:
-
- Allows you to download 3'UTR sequence data directly from biomart for a given subset of genes by manipulating and creating xml strings that are sent to the biomart servers.
- Cleans the 3'UTR data to remove unavailible sequence data and prepares it for direct use within miR-Search.

miR-Search:
-
- Allows you to identify 6mer, 7mer-a1, 7mer-m8 and 8mer targets for a given miRNA.
- Takes a miRNA sequence as an input.
- Searches 3'UTR sequences obtained via pymart for binding sites
- Provides location and type of site