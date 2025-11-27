#!/bin/bash
# Download RefSeq CDS nucleotide sequences for all domains (except bacteria/archaea)
# Uses assembly_summary.txt to find all assemblies and download *_cds_from_genomic.fna.gz

set -e

OUTDIR="/projects/caeg/scratch/kbd606/agp/databases/refseq_cds"
THREADS=16

# Domains to download (skip bacteria/archaea - we have GTDB)
DOMAINS=("fungi" "protozoa" "viral" "invertebrate" "plant" "vertebrate_mammalian" "vertebrate_other")

FTP_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/refseq"

download_domain_cds() {
    local domain=$1
    local outdir="${OUTDIR}/${domain}"

    echo "========================================"
    echo "Processing: ${domain}"
    echo "========================================"

    mkdir -p "${outdir}"
    cd "${outdir}"

    # Download assembly summary
    echo "Downloading assembly_summary.txt..."
    curl -s "${FTP_BASE}/${domain}/assembly_summary.txt" -o assembly_summary.txt

    # Count assemblies
    total=$(grep -v "^#" assembly_summary.txt | wc -l)
    echo "Found ${total} assemblies"

    # Generate CDS URLs from assembly summary
    # Column 20 is ftp_path, we need to append the CDS filename
    echo "Generating CDS URLs..."
    grep -v "^#" assembly_summary.txt | \
        awk -F'\t' '{
            ftp_path=$20;
            if (ftp_path != "na" && ftp_path != "") {
                # Extract assembly name from path
                n=split(ftp_path, parts, "/");
                asm_name=parts[n];
                # Convert https to ftp for aria2
                gsub("https://", "ftp://", ftp_path);
                print ftp_path "/" asm_name "_cds_from_genomic.fna.gz"
            }
        }' > urls_cds.txt

    num_urls=$(wc -l < urls_cds.txt)
    echo "Generated ${num_urls} CDS URLs"

    # Download with aria2
    echo "Downloading CDS files..."
    aria2c -x ${THREADS} -j ${THREADS} \
           --auto-file-renaming=false \
           --allow-overwrite=true \
           --continue=true \
           --max-tries=3 \
           --retry-wait=5 \
           -d "${outdir}/cds" \
           -i urls_cds.txt \
           2>&1 | tee download.log

    # Count downloaded files
    downloaded=$(ls -1 "${outdir}/cds/"*.fna.gz 2>/dev/null | wc -l)
    echo "Downloaded ${downloaded} CDS files for ${domain}"
}

# Create output directory
mkdir -p "${OUTDIR}"

# Process each domain
for domain in "${DOMAINS[@]}"; do
    download_domain_cds "${domain}"
done

echo ""
echo "========================================"
echo "Download complete!"
echo "========================================"
du -sh ${OUTDIR}/*/
