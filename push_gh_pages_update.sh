#!/bin/bash
# Script to push the 404.html update to gh-pages branch
# Run this script from the repository root with appropriate GitHub permissions

set -e

echo "Pushing 404.html update to gh-pages branch..."
echo "Current branch: $(git branch --show-current)"

# Fetch latest changes
git fetch origin gh-pages

# Checkout gh-pages
git checkout gh-pages

# Verify the commit exists
EXPECTED_COMMIT="c7da75613"
CURRENT_COMMIT=$(git rev-parse --short HEAD)
if [[ "$CURRENT_COMMIT" == "$EXPECTED_COMMIT"* ]] || git log --oneline -5 | grep -q "$EXPECTED_COMMIT"; then
    echo "✓ Found the 404.html update commit"
    echo "Commit: $(git log --oneline -1)"
    
    # Push to origin
    echo "Pushing to origin/gh-pages..."
    git push origin gh-pages
    
    echo "✓ Successfully pushed 404.html update to gh-pages!"
    echo ""
    echo "The 404 page should now show version 1.4.0 at:"
    echo "https://pharmaverse.github.io/admiral/404.html"
else
    echo "✗ Could not find the expected commit. Please verify the changes manually."
    git log --oneline -3
    exit 1
fi
