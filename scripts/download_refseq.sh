#!/bin/bash
# Download RefSeq CDS and protein sequences for domain-specific hexamer tables
# Uses aria2 for parallel downloading

set -e

OUTDIR="/projects/caeg/scratch/kbd606/agp/databases/refseq"
THREADS=16

# RefSeq release FTP base
FTP_BASE="ftp://ftp.ncbi.nlm.nih.gov/refseq/release"

# Domains to download
DOMAINS=("fungi" "viral" "archaea" "bacteria" "protozoa")

# File patterns we want:
# - *.rna.fna.gz         - RNA/CDS nucleotide sequences
# - *.protein.faa.gz     - Protein sequences (amino acids)

download_domain() {
    local domain=$1
    local outdir="${OUTDIR}/${domain}"

    echo "=== Downloading ${domain} ==="
    mkdir -p "${outdir}"
    cd "${outdir}"

    # Create URL list for aria2
    # Get directory listing and filter for files we want
    echo "Fetching file list for ${domain}..."

    # Download NT (CDS/RNA) files
    aria2c -x ${THREADS} -j ${THREADS} --auto-file-renaming=false \
           --allow-overwrite=true --continue=true \
           --ftp-pasv=true \
           -i <(lftp -c "open ${FTP_BASE}/${domain}/; cls -1" 2>/dev/null | \
                grep -E '\.(rna\.fna\.gz|cds_from_genomic\.fna\.gz)$' | \
                while read f; do echo "${FTP_BASE}/${domain}/$f"; done) \
           2>&1 | tee "${outdir}/download_nt.log" || true

    # Download AA (protein) files
    aria2c -x ${THREADS} -j ${THREADS} --auto-file-renaming=false \
           --allow-overwrite=true --continue=true \
           --ftp-pasv=true \
           -i <(lftp -c "open ${FTP_BASE}/${domain}/; cls -1" 2>/dev/null | \
                grep -E '\.protein\.faa\.gz$' | \
                while read f; do echo "${FTP_BASE}/${domain}/$f"; done) \
           2>&1 | tee "${outdir}/download_aa.log" || true

    echo "Done with ${domain}"
}

# Check for lftp
if ! command -v lftp &> /dev/null; then
    echo "lftp not found, using alternative method..."
    USE_LFTP=false
else
    USE_LFTP=true
fi

# Download each domain
for domain in "${DOMAINS[@]}"; do
    echo ""
    echo "=========================================="
    echo "Processing: ${domain}"
    echo "=========================================="

    outdir="${OUTDIR}/${domain}"
    mkdir -p "${outdir}"
    cd "${outdir}"

    if [ "$USE_LFTP" = true ]; then
        download_domain "${domain}"
    else
        # Simpler approach: direct aria2 download with FTP directory listing
        # Download NT files
        echo "Downloading ${domain} NT files..."
        aria2c -x ${THREADS} -j ${THREADS} --auto-file-renaming=false \
               --allow-overwrite=true --continue=true \
               "${FTP_BASE}/${domain}/*.rna.fna.gz" \
               2>&1 | tee "download_nt.log" || true

        # Download AA files
        echo "Downloading ${domain} AA files..."
        aria2c -x ${THREADS} -j ${THREADS} --auto-file-renaming=false \
               --allow-overwrite=true --continue=true \
               "${FTP_BASE}/${domain}/*.protein.faa.gz" \
               2>&1 | tee "download_aa.log" || true
    fi
done

echo ""
echo "=== Download complete ==="
echo "Files saved to: ${OUTDIR}"
du -sh ${OUTDIR}/*/
