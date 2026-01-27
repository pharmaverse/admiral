#!/bin/bash
# Script to push the 404.html update branch to GitHub
# Run this script from the repository root with appropriate GitHub permissions

set -e

BRANCH_NAME="copilot/gh-pages-404-update"
EXPECTED_COMMIT="41302eeb4"

echo "Pushing 404.html update branch to GitHub..."

# Checkout the branch
echo "Checking out branch: $BRANCH_NAME"
git checkout "$BRANCH_NAME"

# Verify the commit exists
CURRENT_COMMIT=$(git rev-parse --short HEAD)
if [[ "$CURRENT_COMMIT" == "$EXPECTED_COMMIT"* ]] || git log --oneline -5 | grep -q "$EXPECTED_COMMIT" || true; then
    if git log --oneline -5 | grep -q "$EXPECTED_COMMIT"; then
        echo "✓ Found the 404.html update commit"
        echo "Commit: $(git log --oneline -1)"
        
        # Push to origin
        echo "Pushing $BRANCH_NAME to origin..."
        git push -u origin "$BRANCH_NAME"
        
        echo "✓ Successfully pushed branch $BRANCH_NAME!"
        echo ""
        echo "You can now create a PR to merge this branch into gh-pages."
        echo "Once merged, the 404 page will show version 1.4.0 at:"
        echo "https://pharmaverse.github.io/admiral/404.html"
    else
        echo "✗ Could not find the expected commit. Please verify the changes manually."
        git log --oneline -3
        exit 1
    fi
