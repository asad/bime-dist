#!/usr/bin/env bash
# tools/sign-release.sh - GPG signing helper for BIME release tags.
#
# Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
# Licensed under the Apache License, Version 2.0 - see LICENSE.txt
#
# Opt-in: BIME releases are unsigned by default. Maintainers who want a
# trust-chained release can use this helper to sign (or re-sign) a tag with
# their GPG key, then push it. Consumers verify with `git tag -v <tag>`.
#
# Usage:  bash tools/sign-release.sh v1.1.2
#         bash tools/sign-release.sh             # signs the latest tag
#
# Requires `gpg` and a configured signing key (see `git config user.signingkey`).
# Exits 0 if gpg is missing (with install hints) so this never blocks a CI
# job that simply forgot to install gpg.

set -euo pipefail

if ! command -v gpg >/dev/null 2>&1; then
    echo "gpg not installed - skipping signing."
    echo ""
    echo "  macOS:  brew install gnupg"
    echo "  Debian: sudo apt install gnupg"
    echo "  Fedora: sudo dnf install gnupg2"
    echo ""
    echo "After install, configure git:"
    echo "  git config --global user.signingkey <KEYID>"
    echo "  git config --global commit.gpgsign true"
    exit 0
fi

TAG="${1:-}"
if [ -z "$TAG" ]; then
    TAG=$(git describe --tags --abbrev=0 2>/dev/null || echo "")
    if [ -z "$TAG" ]; then
        echo "no tags in repo and no tag argument given." >&2
        echo "usage: bash tools/sign-release.sh <tag>" >&2
        exit 1
    fi
    echo "no tag given - using latest tag: $TAG"
fi

if ! git rev-parse "$TAG" >/dev/null 2>&1; then
    echo "tag does not exist: $TAG" >&2
    exit 1
fi

# If the tag is already signed and verifies, just confirm.
if git tag -v "$TAG" >/dev/null 2>&1; then
    echo "tag $TAG is already signed and verifies cleanly."
    git tag -v "$TAG"
    exit 0
fi

# Re-sign in place. git refuses to overwrite an existing tag without -f.
echo "signing tag $TAG ..."
COMMIT=$(git rev-list -n 1 "$TAG")
MSG=$(git tag -l --format='%(contents)' "$TAG")
if [ -z "$MSG" ]; then
    MSG="Release $TAG"
fi

git tag -d "$TAG"
git tag -s "$TAG" "$COMMIT" -m "$MSG"

echo ""
echo "verifying signature ..."
git tag -v "$TAG"

echo ""
echo "done. To publish the signed tag:"
echo "  git push --tags --force-with-lease origin $TAG"
echo ""
echo "Consumers verify with:"
echo "  git tag -v $TAG"
