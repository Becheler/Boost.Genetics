# Validation and Benchmarking Plan

## Gold Standard Tool: bcftools

**bcftools** is the industry-standard VCF/BCF manipulation tool from the SAMtools project.

### Installation
```bash
# macOS
brew install bcftools

# Ubuntu/Debian
sudo apt-get install bcftools

# Verify installation
bcftools --version
```

### Why bcftools?
- De facto standard for VCF/BCF processing
- Used by major genomics centers worldwide
- Maintained by original htslib developers
- Comprehensive VCF validation (`bcftools view`, `bcftools query`)

---

## Reference Datasets

### 1. Official VCF Examples (Minimal)
**Source:** VCF v4.3 specification examples
- Location: `tests/data/spec_examples/`
- Small, hand-crafted examples covering edge cases
- Easy to debug, fast to test

### 2. 1000 Genomes Project (Medium)
**Source:** https://www.internationalgenome.org/data
- `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
- Chr22 only (~15 MB compressed, ~500 MB uncompressed)
- 2,504 samples, real-world variants
- Tests: Multi-sample, phased genotypes, INFO fields

### 3. gnomAD (Large - Optional)
**Source:** https://gnomad.broadinstitute.org/
- Subset: gnomAD v3.1 chr22
- Very large, production-scale dataset
- Tests: Performance under load

### 4. Edge Cases Collection
**Custom created for:**
- Structural variants (DEL, INS, DUP, INV, CNV)
- Breakends
- gVCF reference blocks
- Complex genotypes (triploid, mixed ploidy)
- Missing values in various fields
- Large INFO/FORMAT fields

---

## Validation Strategy

### Phase 1: Correctness Tests
Compare our parsing to bcftools on reference datasets:

```bash
# Extract specific fields with bcftools
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n' test.vcf.gz

# Our tool should produce identical output
```

**Test cases:**
1. **Record count**: Same number of variants
2. **Field values**: CHROM, POS, ID, REF, ALT, QUAL, FILTER match exactly
3. **INFO parsing**: Each INFO field matches
4. **Sample genotypes**: GT, DP, GQ, etc. match per sample
5. **Round-trip**: Read → Write → Read produces identical data

### Phase 2: Edge Case Validation
Test specific VCF features:
- SNPs vs Indels vs SVs
- Phased vs unphased genotypes
- Missing values (`.`, `./.`, etc.)
- Multi-allelic sites
- Complex ALT alleles (symbolic, breakends)

### Phase 3: Stress Testing
- Large files (100K+ variants)
- Many samples (1000+)
- Deep nesting in INFO fields
- Memory usage monitoring

---

## Benchmark Plan

### Metrics to Track

**1. Parsing Speed**
- Variants per second
- Records/sec on chr22 (~1M variants)
- Compare to: bcftools, htslib, other C++ libraries

**2. Memory Usage**
- Peak RSS (Resident Set Size)
- Memory per variant
- Header metadata overhead

**3. Write Speed**
- Variants written per second
- Compression ratio (if applicable)

**4. Round-trip Time**
- Read + Write total time

### Benchmark Datasets

| Dataset | Variants | Samples | File Size | Purpose |
|---------|----------|---------|-----------|---------|
| Small | 100 | 1 | ~10 KB | Unit test speed |
| Medium | 10,000 | 10 | ~1 MB | Typical use case |
| Large | 1,000,000 | 100 | ~100 MB | Stress test |
| Chr22 | ~1.1M | 2,504 | ~15 MB gz | Real-world |

### Tools for Comparison

**Baseline (C/C++):**
- bcftools (C + htslib)
- SeqAn3 VCF module
- htslib raw API

**Other languages (informational):**
- cyvcf2 (Python wrapper around htslib)
- VariantWorks (NVIDIA GPU-accelerated)

---

## Implementation Tasks

### 1. Download Reference Data
```bash
mkdir -p tests/data/reference
cd tests/data/reference

# 1000 Genomes chr22
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/\
ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/\
ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
```

### 2. Create Validation Test Suite
`tests/test_validation_reference.cpp`:
- Parse reference VCF
- Compare to bcftools output
- Assert field-by-field equality

### 3. Create Benchmark Suite
`benchmarks/benchmark_parsing.cpp`:
- Time parsing of various datasets
- Memory profiling
- Output: CSV with metrics

### 4. CI Integration
- Run validation tests on every commit
- Performance regression detection
- Compare to previous baseline

---

## Expected Outcomes

**Correctness:**
- ✅ 100% agreement with bcftools on all test datasets
- ✅ Round-trip lossless (read → write → read identical)
- ✅ Pass all VCF v4.3 spec examples

**Performance Targets (relative to bcftools):**
- Parsing: Within 2x of bcftools speed (C library)
- Memory: Similar or better (header-only, no malloc overhead)
- Write: Within 2x of bcftools

**If slower:** Optimization opportunities exist (mmap, SIMD, parallel parsing)

---

## Next Steps

1. ✅ Install bcftools
2. ⬜ Download chr22 test dataset
3. ⬜ Create validation test comparing to bcftools output
4. ⬜ Create benchmark infrastructure
5. ⬜ Run baseline performance measurements
6. ⬜ Optimize hot paths based on profiling
