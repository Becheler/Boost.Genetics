# Comprehensive Benchmark Guide

This document describes the complete benchmark suite for Boost.Genetics, including memory profiling, scalability testing, and writer performance analysis.

## Overview

The benchmark suite consists of five main components:

1. **VCF Parsing Benchmark** (`benchmark_vcf_parsing`) - Sequential parsing performance
2. **Parallel Benchmark** (`benchmark_vcf_parallel`) - Multi-threaded parsing performance
3. **Memory Profiling** (`benchmark_memory`) - Memory usage and arena efficiency
4. **Writer Performance** (`benchmark_writer`) - VCF writing throughput
5. **Scalability Testing** (`benchmark_scalability`) - Performance at different scales

## Building Benchmarks

```bash
cmake -B build -DBUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Running Benchmarks

### Quick Start

Run all benchmarks with a single script:

```bash
./scripts/run_comprehensive_benchmarks.sh tests/data/reference/chr22_subset.vcf
```

### Individual Benchmarks

#### 1. Memory Profiling

Tracks peak memory usage, arena allocation efficiency, and per-record memory overhead:

```bash
./build/benchmark_memory <vcf_file> [num_threads]
```

**Example:**
```bash
./build/benchmark_memory tests/data/reference/chr22_subset.vcf 4
```

**Output metrics:**
- Peak RSS (resident set size)
- Memory overhead per record
- Arena allocator efficiency
- Memory usage for parallel vs sequential

#### 2. Writer Performance

Measures VCF writing throughput and compares with reading:

```bash
# Test with real data
./build/benchmark_writer <input_vcf> [output_vcf]

# Test with synthetic data
./build/benchmark_writer --synthetic <output_vcf> <num_records>
```

**Examples:**
```bash
# Real data
./build/benchmark_writer tests/data/reference/chr22_subset.vcf output.vcf

# Synthetic data - 1 million records
./build/benchmark_writer --synthetic output.vcf 1000000
```

**Output metrics:**
- Read rate (records/sec)
- Write rate (records/sec)
- Write/Read ratio
- Bytes per record
- Write throughput (MB/sec)
- Round-trip performance

#### 3. Scalability Testing

Tests performance consistency across different dataset sizes:

```bash
# Quick test (1K, 10K, 100K records)
./build/benchmark_scalability --quick

# Standard test (10K, 100K, 1M records)
./build/benchmark_scalability

# Full test (includes 10M records - takes longer)
./build/benchmark_scalability --full

# Parallel scalability with existing file
./build/benchmark_scalability --parallel <vcf_file>
```

**Output metrics:**
- Performance at different scales (10K, 100K, 1M, 10M records)
- Scaling consistency (should be O(n))
- Parallel efficiency at different thread counts
- Speedup vs sequential

## Interpreting Results

### Memory Profiling

**Good indicators:**
- Peak RSS scales linearly with file size
- Arena efficiency < 2x overhead
- Bytes per record < 500 bytes

**Example output:**
```
=== Sequential Reader Memory Benchmark ===
Records: 10000
Peak RSS: 45.2 MB
Memory overhead: 12.3 MB
Bytes per record: 1230
```

### Writer Performance

**Good indicators:**
- Write rate > 80% of read rate
- Round-trip performance close to read-only
- Consistent throughput across dataset sizes

**Expected ranges:**
- Sequential write: 3,000-5,000 records/sec
- Write throughput: 5-10 MB/sec

### Scalability

**Good indicators:**
- Rate variation < 20% across dataset sizes (confirms O(n))
- Parallel efficiency > 70% with 4 threads
- Near-linear scaling up to 4 cores

**Example output:**
```
Records     Seq (rec/s)    Par (rec/s)    Speedup
-------     -----------    -----------    -------
10000       4200           16800          4.0x
100000      4150           16600          4.0x
1000000     4180           16720          4.0x

✓ Consistent O(n) scaling (< 20% variation)
```

## Benchmark Suite in CI

The comprehensive benchmarks run automatically in CI:

- **Trigger:** Every push to main, PRs, and manual dispatch
- **Platform:** Ubuntu latest
- **Output:** Benchmark results uploaded as artifacts
- **Reports:** Markdown summary generated

View results in GitHub Actions artifacts after each run.

## Performance Targets

Based on current implementation:

| Metric | Target | Current |
|--------|--------|---------|
| Sequential read | > 4,000 rec/s | ~4,300 rec/s |
| Parallel read (8 threads) | > 15,000 rec/s | ~19,300 rec/s |
| Parallel efficiency (4 threads) | > 70% | ~85% |
| Memory per record | < 500 bytes | ~300 bytes |
| Writer throughput | > 3,000 rec/s | ~3,500 rec/s |
| Scaling consistency | < 20% variation | ~15% variation |

## Profiling for Development

### Using gperftools (optional)

If gperftools is available, memory and CPU profiling are enabled:

```bash
# CPU profiling
CPUPROFILE=prof.out ./build/benchmark_vcf_parsing <file>
google-pprof --text ./build/benchmark_vcf_parsing prof.out

# Memory profiling
HEAPPROFILE=heap.out ./build/benchmark_memory <file>
google-pprof --text ./build/benchmark_memory heap.out
```

### Using perf (Linux)

```bash
perf record -g ./build/benchmark_vcf_parsing <file>
perf report
```

### Using Instruments (macOS)

```bash
instruments -t "Time Profiler" ./build/benchmark_vcf_parsing <file>
```

## Continuous Performance Monitoring

The benchmark suite helps track:

1. **Performance regressions** - Compare with previous runs
2. **Memory efficiency** - Track arena allocator overhead
3. **Scalability** - Ensure O(n) scaling is maintained
4. **Writer performance** - Verify writing doesn't degrade

## Adding New Benchmarks

To add a new benchmark:

1. Create `benchmarks/benchmark_<name>.cpp`
2. Add executable to `CMakeLists.txt`
3. Update `run_comprehensive_benchmarks.sh`
4. Document in this file

## Troubleshooting

**Benchmark fails with "file not found":**
- Ensure test data is downloaded: `./scripts/run_benchmarks.sh --download-data`

**Memory benchmark shows 0 MB:**
- Memory profiling may not be available on your platform
- Check `getCurrentRSS()` implementation for your OS

**Scalability test is slow:**
- Use `--quick` mode for faster testing
- Large datasets (10M records) take several minutes

**Parallel benchmark slower than sequential:**
- Check CPU count with `nproc` or `sysctl -n hw.ncpu`
- Ensure file is large enough (>10K records) for parallelism benefit

## References

- [BENCHMARKING.md](../docs/BENCHMARKING.md) - Quick start guide
- [PERFORMANCE.md](../doc/PERFORMANCE.md) - Performance optimization details
- [GitHub Actions](../.github/workflows/benchmark.yml) - CI benchmark configuration
