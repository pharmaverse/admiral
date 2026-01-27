# 404 Page Update - Completed

## Summary
The 404.html file on the `gh-pages` branch has been successfully updated from version 1.2.0 to 1.4.0.

## Change Details
- **Branch**: `gh-pages`
- **Commit**: `c7da75613`
- **Commit Message**: "Update 404.html to reference version 1.4.0 instead of 1.2.0"

## What Changed
The top-level `404.html` was replaced with the content from `v1.4.0/404.html`, with all URLs adjusted to point to the top-level site instead of the `/v1.4.0/` subdirectory. The version number in the navigation bar now correctly shows **1.4.0** instead of 1.2.0.

## Manual Push Required
The change has been committed locally on the `gh-pages` branch but requires manual pushing with appropriate repository permissions:

```bash
git fetch
git checkout gh-pages  
git push origin gh-pages
```

This will resolve the issue where the 404 page was incorrectly showing version 1.2.0.
