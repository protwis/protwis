
#!/bin/bash
# This script downloads resources from NCBI and generates a lookup file mapping gene symbols and taxonomic IDs to Entrez Gene IDs and species names.

# Function to display usage
usage() {
    echo "Usage: $0 -u uniprot_dir -w working_dir -o outfile"
    echo "  -u uniprot_dir    Directory containing UniProt flat files"
    echo "  -w working_dir    Scratch directory for intermediates and output"
    echo "  -o outfile        Output file path"
    echo "  -t threads        Number of threads to use for parallel processing (default: 10)"
    echo "  -h                Display this help message"
    exit 1
}

threads=10

# Parse command line arguments
while getopts "u:w:o:t:h" opt; do
    case $opt in
        u) uniprot_dir="$OPTARG" ;;
        w) working_dir="$OPTARG" ;;
        o) outfile="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check if required arguments are provided
if [ -z "$uniprot_dir" ] || [ -z "$working_dir" ] || [ -z "$outfile" ]; then
    echo "Error: Missing required arguments."
    usage
fi

#download the gene2accession file from NCBI, which contains mappings between gene symbols and Entrez Gene IDs, as well as taxonomic information
if [ -e "${working_dir}/gene2accession.gz" ]
then
    echo "Found gene2accession.gz in working directory, skipping download. If you want to redownload, please remove the existing file and rerun the script."
else
    echo "Retrieving gene2accession.gz from NCBI FTP server..."
    curl ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz -o ${working_dir}/gene2accession.gz
fi

#download the taxonomy data dump files which various files with taxonomic information, including the mapping between taxonomic IDs and species names
if [ -e "${working_dir}/taxdump.tar.gz" ]
then
    echo "Found taxdump.tar.gz in working directory, skipping download. If you want to redownload, please remove the existing file and rerun the script."
else
    echo "Retrieving taxdump.tar.gz from NCBI FTP server..."
    curl ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -o ${working_dir}/taxdump.tar.gz
fi

#extract
echo "Extracting taxdump.tar.gz..."
cd ${working_dir} && tar -xf taxdump.tar.gz

mkdir -p ${working_dir}/batches
mkdir -p ${working_dir}/batches/processed

#Generate batches of 5000 files for parallel processing.
find ${uniprot_dir} -iname "*.txt" | sort | split -l 5000 - ${working_dir}/batches/uniprot_batch_ --additional-suffix=.txt

#define function to extract taxonomic ID and gene symbol lines from uniprot files, and export it for use with parallel
extract_id_lines() {
        infile=$1
        outfile=$(dirname ${infile})/processed/$(basename ${infile} .txt)_idlines.txt
        cat ${infile}|
        while read -r f;
        do
            grep -E "^OX\s+NCBI_TaxID|^GN\s+Name|^GN\s+Synonyms|DR\s+GeneID;" ${f};
        done > ${outfile}
}
export -f extract_id_lines

#Use parallel to process batches
parallel --jobs ${threads} extract_id_lines ::: ${working_dir}/batches/uniprot_batch_*.txt

#Merge batch outputs
cat ${working_dir}/batches/processed/*_idlines.txt > ${working_dir}/uniprot_taxid_genesymbol_intermediate.txt

#Use awk to create three intermediate files: one with taxonomic IDs, one with unique gene symbols (including synonyms), and one with Entrez Gene IDs from the uniprot files.
awk \
    -v FILEa=${working_dir}/unique_taxid_list.txt \
    -v FILEb=${working_dir}/unique_genesymbol_list.txt \
    -v FILEc=${working_dir}/unique_entrezgeneid_list.txt \
    'BEGIN{
        FS="\t"; OFS="\t"
    }

    /^OX/{
        match($0, "^OX[[:space:]]+NCBI_TaxID=([0-9]+) ?[{;].*$", match_array);
        entrez_taxid = match_array[1];
        if (entrez_taxid != "") { print entrez_taxid >> FILEa }
        if (entrez_taxid == "") { print ($0, "- taxid not parsed") }
    }

    /^DR   GeneID;/{
        match($0, "^DR[[:space:]]+GeneID;[[:space:]]+([0-9]+)[[:space:]]?;.*", match_array);
        entrez_geneid = match_array[1];
        if (entrez_geneid != "") { print entrez_geneid >> FILEc }
        if (entrez_geneid == "") { print ($0, "- entrez_geneid not parsed") }
    }

    /^GN[[:space:]]{3}Name=/{
        match($0, "^GN[[:space:]]+Name=([^[:space:]{;]+)[[:space:]]?[{;].*$", match_array)
        gene_symbol = match_array[1];
        if (gene_symbol != "") { print gene_symbol >> FILEb }
        if (gene_symbol == "") { print ($0, "- gene_symbol not parsed") }
    }

    /^GN.+Synonyms=/{
        match($0, "^GN.+Synonyms=([^[:space:]{;]+)[^;]*[[:space:]]?;.*", match_array)
        synonyms = match_array[1];
        split(synonyms, syn_array, ", ");
        for (i in syn_array) {
                if (syn_array[i] != "") { print syn_array[i] >> FILEb }
                if (syn_array[i] == "") { print ($0, "- gene_symbol_synonym not parsed") }
            }
    }' ${working_dir}/uniprot_taxid_genesymbol_intermediate.txt

# Sort and remove duplicates from the three intermediate files to create final lists of unique taxonomic IDs, gene symbols,
# and Entrez Gene IDs that we will use to filter the gene2accession file.
sort -n ${working_dir}/unique_taxid_list.txt | uniq > temp.txt && mv temp.txt ${working_dir}/unique_taxid_list.txt
sort ${working_dir}/unique_genesymbol_list.txt | tr -d ',' | uniq > temp.txt && mv temp.txt ${working_dir}/unique_genesymbol_list.txt
sort -n ${working_dir}/unique_entrezgeneid_list.txt | uniq > temp.txt && mv temp.txt ${working_dir}/unique_entrezgeneid_list.txt

#Account for incorrect human gene name in uniprot file.
echo "NPY6RP" >> ${working_dir}/unique_genesymbol_list.txt #Uniprot file for Q99463 currently uses outdated NPY6R gene symbol instead of corrrect NPY6RP

#Extract the mapping between taxonomic IDs and species names from the taxonomy dump file
echo "Extracting scientific names from taxonomy dump..."
grep scientific names.dmp |
	sed "s/\t//g" |
	cut -f1,2 -d'|' |
	tr "|" "\t" > ${working_dir}/taxid_speciesname_mapping.txt

# Compute the number of lines in the gene2accession file and determine how many lines to process per thread based on the specified number of threads.
# We add a buffer of 10,000 lines to ensure that we don't have a small final batch or miss any records.
echo "Computing batch sizes for parallel processing..."
entrez_lines=$(zcat ${working_dir}/gene2accession.gz | wc -l)
lines_per_thread=$(( ((entrez_lines) / threads) + 10000 ))

# Split the gene2accession file into batches for parallel processing
echo "Splitting gene2accession file into batches of ${lines_per_thread} lines for parallel processing..."
zcat ${working_dir}/gene2accession.gz | split -l ${lines_per_thread} -d --additional-suffix=.txt - ${working_dir}/batches/gene2accession_batch_

#Embedded python script to perform the batched filtering of the gene2accession file.
# The script will output a line when it finds a record that matches either a gene symbol and taxonomic ID combination from the uniprot files,
# or an Entrez Gene ID from the uniprot files.
python3 <<-EOF
from concurrent.futures import ThreadPoolExecutor
import glob

print('Loading tax ids, gene symbols and entrez gene ids of interest ...')
with open("${working_dir}/unique_taxid_list.txt") as tax_fh:
    tax_ids = set([line.strip() for line in tax_fh])

with open("${working_dir}/unique_genesymbol_list.txt") as symbol_fh:
    gene_symbols = set([line.strip().strip() for line in symbol_fh])

with open("${working_dir}/unique_entrezgeneid_list.txt") as geneid_fh:
    entrez_ids = set([line.strip() for line in geneid_fh])

print('Loading tax id to specimen lookup ...')
tax_id_dict = dict()
with open("${working_dir}/taxid_speciesname_mapping.txt") as tax_name_fh:
    for line in tax_name_fh:
        tax_id, species_name = line.strip().split('\t')
        tax_id_dict[tax_id] = species_name

print('Filtering Entrez Ids ...')

def process_file(file_path, gene_symbols, tax_ids, entrez_ids, tax_id_dict):
    outfile_fh = open(f"{file_path}.filtered", "w")
    processed_record_counter = 0
    written_record_counter = 0
    with open(file_path, "r") as entrez_fh:
        for line in entrez_fh:
            line_split = line.strip().split('\t')
            if len(line_split) < 16:
                continue
            tax_id = line_split[0]
            gene_id = line_split[1]
            symbol = line_split[15]

            if (symbol in gene_symbols and tax_id in tax_ids) or gene_id in entrez_ids:
                species = tax_id_dict[tax_id]
                outfile_fh.write(f'{symbol}\t{tax_id}\t{gene_id}\t{species}\n')
                written_record_counter += 1

            processed_record_counter += 1
            if processed_record_counter % 1000000 == 0:
                print(f'Processed:{processed_record_counter}, wrote:{written_record_counter} from {file_path}...')
    outfile_fh.close()
    print(f'Finished processing {file_path}. Total processed:{processed_record_counter}, total written:{written_record_counter}.')

files_to_process = glob.glob('${working_dir}/batches/gene2accession_batch_*.txt')
with ThreadPoolExecutor(max_workers=${threads}) as executor:
    for file_path in files_to_process:
        executor.submit(process_file, file_path, gene_symbols, tax_ids, entrez_ids, tax_id_dict)
EOF

echo "Combining batch outputs, sorting and removing duplicates to create final output file..."

#Combine batched output, sort, and remove duplicates
echo -e "gene_symbol\ttaxon_id\tentrez_gene_id\tspecies_name" > ${outfile}
cat ${working_dir}/batches/gene2accession_batch_*.txt.filtered | sort -k3n -k2n | uniq >> ${outfile}

#Account for incorrect human gene name in uniprot file.
echo "NPY6R#9606#4888#Homo sapiens" | tr '#' '\t' >> ${outfile} #Uniprot file for Q99463 currently uses outdated NPY6R gene symbol instead of corrrect NPY6RP

echo "Cleaning up ..."
if [ -e "${outfile}" ]
then
#cleanup batch files
    rm -fr ${working_dir}/batches

    rm ${working_dir}/*.dmp \
    ${working_dir}/*.prt \
    ${working_dir}/readme.txt \
    ${working_dir}/taxid_speciesname_mapping.txt \
    ${working_dir}/unique_genesymbol_list.txt \
    ${working_dir}/unique_taxid_list.txt \
    ${working_dir}/unique_entrezgeneid_list.txt \
    ${working_dir}/uniprot_taxid_genesymbol_intermediate.txt
fi