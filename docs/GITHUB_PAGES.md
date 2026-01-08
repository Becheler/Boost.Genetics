# GitHub Pages Setup

This document explains how to enable GitHub Pages for the Boost.Genetics documentation.

## Enabling GitHub Pages

1. Go to your repository on GitHub
2. Click **Settings** → **Pages** (in the left sidebar)
3. Under **Build and deployment**:
   - **Source**: Select "GitHub Actions"
4. Click **Save**

That's it! The documentation will be automatically deployed on every push to `main`.

## What Gets Deployed

The GitHub Pages site includes:

- **Home page** (`index.html`): Project overview with live benchmark results
- **README** (`README.html`): Full project documentation
- **Benchmarking Guide** (`BENCHMARKING.html`): How to run benchmarks locally
- **Latest Benchmark Results** (`benchmark_results.html`): Auto-updated performance metrics

## Workflows

### `benchmark.yml`
- **Triggers**: Push to main/feature branches, PRs, manual dispatch
- **What it does**:
  - Downloads 1000 Genomes chr22 test dataset
  - Builds benchmark executables
  - Runs sequential and parallel benchmarks
  - Compares with bcftools
  - Posts results as PR comment (on PRs)
  - Uploads results as artifacts

### `pages.yml`
- **Triggers**: Push to main, manual dispatch
- **What it does**:
  - Runs all benchmarks
  - Generates HTML documentation from markdown
  - Creates a static website with live benchmark data
  - Deploys to GitHub Pages

## Viewing Results

Once enabled, your documentation will be available at:
```
https://<username>.github.io/Boost.Genetics/
```

For example:
```
https://arnaudbecheler.github.io/Boost.Genetics/
```

## Badge URLs

Add these to your README for status badges:

```markdown
[![CI](https://github.com/arnaudbecheler/Boost.Genetics/actions/workflows/ci.yml/badge.svg)](https://github.com/arnaudbecheler/Boost.Genetics/actions/workflows/ci.yml)
[![Benchmark](https://github.com/arnaudbecheler/Boost.Genetics/actions/workflows/benchmark.yml/badge.svg)](https://github.com/arnaudbecheler/Boost.Genetics/actions/workflows/benchmark.yml)
[![Documentation](https://github.com/arnaudbecheler/Boost.Genetics/actions/workflows/pages.yml/badge.svg)](https://arnaudbecheler.github.io/Boost.Genetics/)
```

## Local Preview

To preview the documentation locally:

```bash
# Run benchmarks and generate docs
cd Boost.Genetics
mkdir -p _site

# Build benchmarks
cmake -B build -DBUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# Run benchmarks (requires test data)
./scripts/run_benchmarks.sh

# Open in browser
open _site/index.html  # macOS
xdg-open _site/index.html  # Linux
```

## Customization

Edit these files to customize the documentation:

- `.github/workflows/pages.yml` - Deployment workflow
- HTML template in pages.yml (search for `cat > _site/index.html`)
- `docs/BENCHMARKING.md` - Benchmarking guide content
- `README.md` - Main documentation

## Troubleshooting

### Pages not deploying
- Check Actions tab for workflow errors
- Verify Pages is enabled in Settings
- Ensure `pages.yml` workflow completed successfully

### Benchmark failures
- Check if test data downloaded correctly
- Verify bcftools is installed (optional but recommended)
- Check build logs for compilation errors

### Outdated benchmark results
- GitHub Pages caches aggressively
- Add `?nocache=<timestamp>` to URL to bypass
- Wait a few minutes for CDN propagation
