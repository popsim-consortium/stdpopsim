#!/usr/bin/env bash
#
# Upload a tarball to the stdpopsim S3 bucket via GitHub Actions.
#
# Usage:
#   ./maintenance/upload_to_s3.sh <local-file> <species-id> <genetic_map|annotation> <resource-id>
#
# Examples:
#   ./maintenance/upload_to_s3.sh data/HapMapII_GRCh38.tar.gz HomSap genetic_map HapMapII_GRCh38
#   ./maintenance/upload_to_s3.sh data/ensembl_havana_104_exons_v1.tar.gz HomSap annotation ensembl_havana_104_exons
#
# The script looks up the resource in the stdpopsim catalog to find the
# expected S3 URL and SHA256, verifies the local file matches, then uploads
# via a GitHub Actions workflow.
#
# Prerequisites:
#   - GitHub CLI (gh) installed and authenticated
#   - stdpopsim installed in current Python environment (pip install -e .)
#   - Push access to the stdpopsim repository

set -euo pipefail

REPO="popsim-consortium/stdpopsim"
WORKFLOW="upload-to-s3.yml"

die() { echo "ERROR: $*" >&2; exit 1; }

# --- Validate arguments ---
if [ $# -ne 4 ]; then
    echo "Usage: $0 <local-file> <species-id> <genetic_map|annotation> <resource-id>" >&2
    echo "" >&2
    echo "Examples:" >&2
    echo "  $0 data/HapMapII_GRCh38.tar.gz HomSap genetic_map HapMapII_GRCh38" >&2
    echo "  $0 data/exons_v1.tar.gz HomSap annotation ensembl_havana_104_exons" >&2
    exit 1
fi

LOCAL_FILE="$1"
SPECIES_ID="$2"
RESOURCE_TYPE="$3"
RESOURCE_ID="$4"

[ -f "$LOCAL_FILE" ] || die "File not found: $LOCAL_FILE"

case "$LOCAL_FILE" in
    *.tar.gz|*.tgz) ;;
    *) die "File must be a .tar.gz or .tgz archive: $LOCAL_FILE" ;;
esac

case "$RESOURCE_TYPE" in
    genetic_map|annotation) ;;
    *) die "Resource type must be 'genetic_map' or 'annotation', got: $RESOURCE_TYPE" ;;
esac

# --- Check prerequisites ---
command -v gh >/dev/null 2>&1 || die "GitHub CLI (gh) is not installed. See https://cli.github.com/"
gh auth status >/dev/null 2>&1 || die "GitHub CLI is not authenticated. Run: gh auth login"
python3 -c "import stdpopsim" 2>/dev/null || die "stdpopsim is not installed in the current Python environment"

# --- Compute SHA256 ---
if command -v sha256sum >/dev/null 2>&1; then
    SHA256=$(sha256sum "$LOCAL_FILE" | awk '{print $1}')
elif command -v shasum >/dev/null 2>&1; then
    SHA256=$(shasum -a 256 "$LOCAL_FILE" | awk '{print $1}')
else
    die "Neither sha256sum nor shasum found"
fi

# --- Validate against stdpopsim catalog ---
echo "Validating against stdpopsim catalog..."
VALIDATION=$(python3 << PYEOF
import stdpopsim
import sys
import json
import re

species_id = "${SPECIES_ID}"
resource_type = "${RESOURCE_TYPE}"
resource_id = "${RESOURCE_ID}"
local_sha256 = "${SHA256}"

try:
    species = stdpopsim.get_species(species_id)
except (ValueError, KeyError) as e:
    print(json.dumps({"error": f"Species '{species_id}' not found in catalog. "
                       f"Available: {[s.id for s in stdpopsim.all_species()]}"}))
    sys.exit(0)

if resource_type == "genetic_map":
    resources = {gm.id: gm for gm in species.genetic_maps}
    if resource_id not in resources:
        avail = list(resources.keys())
        print(json.dumps({"error": f"Genetic map '{resource_id}' not found for "
                           f"{species_id}. Available: {avail}"}))
        sys.exit(0)
    resource = resources[resource_id]
    s3_url = resource.url
    expected_sha256 = resource.sha256

elif resource_type == "annotation":
    resources = {a.id: a for a in species.annotations}
    if resource_id not in resources:
        avail = list(resources.keys())
        print(json.dumps({"error": f"Annotation '{resource_id}' not found for "
                           f"{species_id}. Available: {avail}"}))
        sys.exit(0)
    resource = resources[resource_id]
    s3_url = resource.intervals_url
    expected_sha256 = resource.intervals_sha256

# Extract S3 path from URL (handles both s3-us-west-2 and s3.us-west-2 styles)
m = re.match(r"https://stdpopsim\.s3[.-]us-west-2\.amazonaws\.com/(.*)", s3_url)
if not m:
    print(json.dumps({"error": f"Could not parse S3 URL from catalog: {s3_url}"}))
    sys.exit(0)

s3_dest = m.group(1)

if local_sha256 != expected_sha256:
    print(json.dumps({"error": f"SHA256 mismatch!\n"
                       f"  Local file:  {local_sha256}\n"
                       f"  Catalog expects: {expected_sha256}"}))
    sys.exit(0)

print(json.dumps({
    "s3_url": s3_url,
    "s3_dest": s3_dest,
    "sha256": expected_sha256,
}))
PYEOF
) || die "Python validation script failed"

# Check for validation error
ERROR=$(echo "$VALIDATION" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('error',''))")
if [ -n "$ERROR" ]; then
    die "$ERROR"
fi

# Extract validated fields
S3_DEST=$(echo "$VALIDATION" | python3 -c "import sys,json; print(json.load(sys.stdin)['s3_dest'])")
S3_URL=$(echo "$VALIDATION" | python3 -c "import sys,json; print(json.load(sys.stdin)['s3_url'])")

echo ""
echo "Catalog validation passed!"
echo "  Species:  $SPECIES_ID"
echo "  Type:     $RESOURCE_TYPE"
echo "  Resource: $RESOURCE_ID"
echo "  File:     $LOCAL_FILE"
echo "  S3 URL:   $S3_URL"
echo "  SHA256:   $SHA256"
echo ""

# --- Create draft release as file transport ---
TAG="s3-upload-$(date +%Y%m%d%H%M%S)-$RANDOM"
FILENAME=$(basename "$LOCAL_FILE")

echo "Creating draft release ($TAG) to transport the file..."
gh release create "$TAG" \
    --repo "$REPO" \
    --draft \
    --title "S3 upload: $SPECIES_ID/$RESOURCE_TYPE/$RESOURCE_ID" \
    --notes "Temporary draft release for S3 upload. Will be cleaned up automatically." \
    "$LOCAL_FILE"

echo "Draft release created."
echo ""

# --- Trigger the workflow ---
echo "Triggering upload workflow..."
gh workflow run "$WORKFLOW" \
    --repo "$REPO" \
    --field "release_tag=$TAG" \
    --field "s3_destination=$S3_DEST" \
    --field "expected_sha256=$SHA256" \
    --field "species_id=$SPECIES_ID" \
    --field "resource_type=$RESOURCE_TYPE" \
    --field "resource_id=$RESOURCE_ID"

echo ""
echo "Workflow triggered successfully."
echo ""
echo "Monitor the upload with:"
echo "  gh run list --repo $REPO --workflow $WORKFLOW --limit 3"
echo ""
echo "Or watch the latest run:"
echo "  gh run watch --repo $REPO \$(gh run list --repo $REPO --workflow $WORKFLOW --limit 1 --json databaseId --jq '.[0].databaseId')"
echo ""
echo "Final S3 URL: $S3_URL"
