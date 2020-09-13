# Bioinformatics-Tools
These files were created as I worked with the Leung Lab Database

cleaning_Domain_text.R is a file that takes a text file downloaded from the Uniprot website. Information about domians for each protien in the file is extracted 
and the Names, details and amino acid positions of each domain is reported. This code was developed to add domian information to a protien graphic, so the names 
and locations of domians were essential.

Sequence_Check.R onfirms single amino acid modifications that I use to verify information in the Leung Lab database. Three files are required: 
text file listing the Uniprot IDs for proteins of interest, a text file listing the respective modifications for each protein of interest, 
and a FASTA file containing information on each protein of interest. Proteins for which the reported modifications are incorrect are returned.

findUniprotLength.R uses a text file from the Uniprot website, and reports the Uniprot ID of proteins listed and their respective amino acid lengths.
