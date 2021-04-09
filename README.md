# miR-Search

miR-Search is a programme developed to identify all 6mer, 7mer and 8mer targets of a given micro-RNA as well as provoiding an interface to download and clean 3'UTR data from the ensembl database to be used in the programme

This project consists of two modules pymart and miR-Search:

pymart:
-
- Allows you to download 3'UTR sequence data directly from biomart for a given subset of genes by manipulating and creating xml strings that are sent to the biomart servers.
- Cleans the 3'UTR data to remove unavailible sequence data and prepares it for direct use within miR-Search.

    
    Usage python miR-Search -d -i "genelist" - o "output" -a

Note: The genelist should be a csv file with the first column the ensembl gene IDs

template.xml contains the xml request string. This currently obtains 3'UTR sequences (specified via: `<Attribute name = "3utr" />`).

### Arguments
- --download, -d: Runs the pymart downloader. Requires --input.
- --input, -i: Specifies the location of gene list in csv format.
- --output, -o: Specifies the location of the FASTA output
- --check, -c: Will run a comparison check between the output gene sequence and the input gene list.
- --auto, -a: Runs pymart in automode, performing the comparison check and file tidy (same as running -c -s)
- --threads, -t: Specifies the number of threads pymart will use when downloading. Default = 50
- --scrub, -s: Will clean the specified FASTA file. Recomended to be run after downloading.
- --url: 

miR-Search:
-
- Allows you to identify 6mer, 7mer-a1, 7mer-m8 and 8mer targets for a given miRNA.
- Takes a miRNA sequence as an input.
- Searches 3'UTR sequences obtained via pymart for binding sites
- Provides location and type of site