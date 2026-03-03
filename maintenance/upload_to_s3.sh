#!/usr/bin/env bash
#
# Upload a tarball to the stdpopsim S3 bucket via GitHub Actions.
#
# Usage:
#   ./maintenance/upload_to_s3.sh <local-file> <s3-destination>
#
# Example:
#   ./maintenance/upload_to_s3.sh annotations/HomSap/exons_v2.tar.gz annotations/HomSap/exons_v2.tar.gz
#
# The s3-destination is the path within the stdpopsim bucket, and must start
# with "genetic_maps/" or "annotations/".
#
# Prerequisites:
#   - GitHub CLI (gh) installed and authenticated
#   - Push access to the stdpopsim repository

set -euo pipefail

REPO="popsim-consortium/stdpopsim"
WORKFLOW="upload-to-s3.yml"

die() { echo "ERROR: $*" >&2; exit 1; }

# --- Validate arguments ---
if [ $# -ne 2 ]; then
    echo "Usage: $0 <local-file> <s3-destination>" >&2
    echo "Example: $0 annotations/HomSap/exons_v2.tar.gz annotations/HomSap/exons_v2.tar.gz" >&2
    exit 1
fi

LOCAL_FILE="$1"
S3_DEST="$2"

[ -f "$LOCAL_FILE" ] || die "File not found: $LOCAL_FILE"

case "$LOCAL_FILE" in
    *.tar.gz) ;;
    *) die "File must be a .tar.gz archive: $LOCAL_FILE" ;;
esac

case "$S3_DEST" in
    genetic_maps/*|annotations/*) ;;
    *) die "S3 destination must start with 'genetic_maps/' or 'annotations/': $S3_DEST" ;;
esac

case "$S3_DEST" in
    *.tar.gz) ;;
    *) die "S3 destination must end with .tar.gz: $S3_DEST" ;;
esac

# --- Check gh CLI ---
command -v gh >/dev/null 2>&1 || die "GitHub CLI (gh) is not installed. See https://cli.github.com/"
gh auth status >/dev/null 2>&1 || die "GitHub CLI is not authenticated. Run: gh auth login"

# --- Compute SHA256 ---
if command -v sha256sum >/dev/null 2>&1; then
    SHA256=$(sha256sum "$LOCAL_FILE" | awk '{print $1}')
elif command -v shasum >/dev/null 2>&1; then
    SHA256=$(shasum -a 256 "$LOCAL_FILE" | awk '{print $1}')
else
    die "Neither sha256sum nor shasum found"
fi

echo "File:        $LOCAL_FILE"
echo "Destination: s3://stdpopsim/$S3_DEST"
echo "SHA256:      $SHA256"
echo ""

# --- Create draft release as file transport ---
TAG="s3-upload-$(date +%Y%m%d%H%M%S)-$RANDOM"
FILENAME=$(basename "$LOCAL_FILE")

echo "Creating draft release ($TAG) to transport the file..."
gh release create "$TAG" \
    --repo "$REPO" \
    --draft \
    --title "S3 upload: $FILENAME" \
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
    --field "expected_sha256=$SHA256"

echo ""
echo "Workflow triggered successfully."
echo ""
echo "Monitor the upload with:"
echo "  gh run list --repo $REPO --workflow $WORKFLOW --limit 3"
echo ""
echo "Or watch the latest run:"
echo "  gh run watch --repo $REPO \$(gh run list --repo $REPO --workflow $WORKFLOW --limit 1 --json databaseId --jq '.[0].databaseId')"
echo ""
echo "Final S3 URL: https://stdpopsim.s3-us-west-2.amazonaws.com/$S3_DEST"
echo ""
echo "SHA256 for catalog code: $SHA256"
