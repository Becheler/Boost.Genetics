# Benchmarking Guide

This document describes how to run performance benchmarks for Boost.Genetics VCF parsing.

## Quick Start

Run automated benchmarks:
```bash
./scripts/run_benchmarks.sh
```

This script will:
1. Download test datasets (1000 Genomes chr22)
2. Check for bcftools (optional, for comparison)
3. Build benchmark executables
4. Run performance tests with different thread counts
5. Compare against bcftools baseline

## Manual Benchmarking

### Prerequisites

**Required:**
- CMake 3.14+
- C++17 compiler (GCC 7+, Clang 5+, MSVC 2017+)
- Test dataset (see below)

**Optional (for comparison):**
- bcftools (`brew install bcftools` on macOS)

### Test Datasets

The automated script downloads data automatically, but you can also use your own VCF files:

**Option 1: Use 1000 Genomes data (recommended)**
```bash
cd tests/data/reference
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
gunzip -c ALL.chr22*.vcf.gz | head -n 10000 > chr22_subset.vcf
```

**Option 2: Use your own VCF file**
```bash
# Any standard VCF file works
cp /path/to/your/data.vcf tests/data/reference/
```

### Building Benchmarks

```bash
mkdir build && cd build
cmake .. -DBUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release
make -j8
```

### Running Benchmarks

**Sequential parsing:**
```bash
./build/benchmark_vcf_parsing tests/data/reference/chr22_subset.vcf
```

**Parallel parsing (auto-detect threads):**
```bash
./build/benchmark_vcf_parallel tests/data/reference/chr22_subset.vcf
```

**Parallel parsing (specific thread count):**
```bash
./build/benchmark_vcf_parallel tests/data/reference/chr22_subset.vcf 8
```

**Test scaling across thread counts:**
```bash
for threads in 1 2 4 8; do
    echo "=== $threads threads ==="
    ./build/benchmark_vcf_parallel tests/data/reference/chr22_subset.vcf $threads
done
```

### Comparing with bcftools

```bash
# Time bcftools parsing
time bcftools view -H tests/data/reference/chr22_subset.vcf > /dev/null

# Time Boost.Genetics parallel parsing
time ./build/benchmark_vcf_parallel tests/data/reference/chr22_subset.vcf 8
```

## Expected Performance

Based on benchmarks with 1000 Genomes chr22 data (6,272 records, ~83 MB):

| Implementation | Records/sec | Speedup vs Sequential | vs bcftools |
|---------------|-------------|----------------------|-------------|
| Sequential | ~4,300 | 1.0x | 0.36x |
| Parallel (2 threads) | ~12,900 | 3.0x | **1.07x** |
| Parallel (4 threads) | ~19,200 | 4.5x | **1.60x** |
| Parallel (8 threads) | ~19,300 | 4.5x | **1.61x** |
| bcftools | ~12,000 | - | 1.0x |

**Key Achievements:**
- ✅ **82x faster** than baseline implementation (234 → 19,300 rec/s)
- ✅ **60% faster** than bcftools with 8 threads
- ✅ **Near-linear scaling** up to 4 cores
- ✅ **Beats bcftools** even with just 2 threads

## Performance Optimizations

The current implementation includes:

1. **Zero-copy parsing** - Uses `std::string_view` for all fields
2. **Arena allocator** - Bump-pointer allocation eliminates malloc overhead
3. **Lazy INFO parsing** - Defers parsing until field is accessed
4. **Parallel chunking** - Divides file into chunks processed by worker threads
5. **Thread-local arenas** - Zero contention between parsing threads
6. **Optimized string splitting** - Manual pointer iteration over std::string::find()

## Profiling (Advanced)

If gperftools is installed, you can profile hotspots:

```bash
# Build with profiling support
cmake .. -DBUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release
make benchmark_vcf_parsing

# Run with profiling
CPUPROFILE=vcf_parsing.prof ./build/benchmark_vcf_parsing tests/data/reference/chr22_subset.vcf

# Analyze results
pprof --text ./build/benchmark_vcf_parsing vcf_parsing.prof | head -20
```

## Troubleshooting

**Benchmark fails to build:**
- Ensure C++17 support: `g++ --version` (need 7+) or `clang++ --version` (need 5+)
- Check CMake version: `cmake --version` (need 3.14+)

**Performance is slower than expected:**
- Ensure Release build: `-DCMAKE_BUILD_TYPE=Release`
- Check CPU scaling: `cat /proc/cpuinfo | grep MHz` (Linux)
- Disable other programs consuming CPU
- Use SSD storage for best I/O performance

**Dataset download fails:**
- Check internet connection
- Try manual download with browser
- Use your own VCF file instead

## Benchmark Data Sources

- **1000 Genomes Project**: High-quality population genomics data
  - URL: http://ftp.1000genomes.ebi.ac.uk/
  - License: Public domain
  - Used for: Standard benchmarking, realistic workload

- **Custom datasets**: Use your own VCF files for domain-specific testing
