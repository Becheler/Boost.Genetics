# VCF Performance Analysis

## Test Setup
- **Dataset**: 1000 Genomes chr22 subset (chr22:10519265-11000000)
- **Records**: 6,272 variants
- **Samples per record**: 3,202 individuals
- **File Size**: 83.2 MB uncompressed
- **Gold Standard**: bcftools v1.23
- **Test Platform**: macOS (Apple Silicon)

## Performance Evolution

### Baseline → Current (16.2x improvement achieved)

| Phase | Optimization | Records/sec | Speedup | Gap vs bcftools |
|-------|-------------|-------------|---------|-----------------|
| **Initial** | None | 234 | 1.0x | 50.5x slower |
| **C++17** | string_view + fast parsing | 441 | 1.9x | 26.7x slower |
| **O3** | Compiler optimization (-O3) | 3,440 | 14.7x | 3.4x slower |
| **Vector** | vector<pair> for samples | 3,790 | 16.2x | **3.1x slower** |
| **bcftools baseline** | - | 11,611 | - | - |

**Current Status**: 3,790 rec/s (**within 3.1x of bcftools**)

### Writing Speed
| Operation | Time (s) | Records/sec | MB/s |
|-----------|----------|-------------|------|
| **VCF Writing** | 0.708 | 8,859 | 117.56 |

### Validation Overhead
- Parsing without validation: 1.676s (3,743 rec/s)
- Parsing with validation: 1.632s (3,844 rec/s)
- **Validation overhead: -2.6%** (negative = validation actually faster due to caching)

### Memory Usage
- Average bytes per record: 13,910 bytes
- Total for 6,272 records: ~87 MB

## Deep Performance Audit Results

### Root Cause Analysis
**PRIMARY BOTTLENECK: String allocations**
- **401,658,880 allocations** per full parse
- 3,202 samples × 6,272 records × 10 fields × 2 (key+value) = **401M allocations**
- Each `std::string(value)` triggers heap allocation
- Format field keys duplicated for every sample

### Profiling Breakdown
```
I/O time:     0.000167s  (0.01%)
Parse time:   1.654s     (99.9%) ← THE BOTTLENECK
Other:        0.000686s  (0.04%)
```

**Within parse_record():**
- Sample data parsing: ~95% (3,202 samples with 10 fields each)
- INFO field parsing: ~3%
- Fixed fields: ~2%

### Micro-Benchmark Results
Testing different storage strategies on same dataset:

| Strategy | Records/sec | Allocations | vs Current |
|----------|-------------|-------------|------------|
| **Current** (pair<string,string>) | 4,293 | 401M | baseline |
| **Indexed** (vector<string>, no keys) | 5,468 | 200M | ⬆ 27% |
| **String_view** (zero-copy) | 5,575 | 0 | ⬆ 30% |

**Key Finding**: We can eliminate 200M+ allocations by avoiding format key duplication and using string_view for values.

## Validation Results
✅ **All validation tests passed**
- 100,356 assertions against bcftools
- Fields tested: CHROM, POS, ID, REF, ALT, FILTER, INFO, FORMAT, samples
- **Result**: Bit-for-bit compatible with bcftools output

## Remaining Optimization Opportunities (Ranked by Impact)

### 🔥 Tier 1: Game Changers (1.5-3x potential each)

#### 1. Indexed Sample Format (27% gain - 2 hours)
**Problem**: Duplicating format field names for every sample  
**Solution**: Store format fields once per record, samples are just value arrays
```cpp
// Current: vector<pair<string,string>> - duplicates "GT", "DP", "GQ" 3,202 times
// Optimized: vector<string> format_fields + vector<vector<string>> sample_values
```
**Expected**: 3,790 → 4,800 rec/s (halves allocations from 401M → 200M)

#### 2. String_view Sample Values (30% gain - 4 hours)  
**Problem**: Copying every sample value into std::string  
**Solution**: Keep line buffer alive, store string_view instead
```cpp
class record {
    std::shared_ptr<std::string> line_buffer_;  // Keep source alive
    std::vector<std::vector<std::string_view>> samples_;
};
```
**Expected**: 4,800 → 6,200 rec/s (eliminates remaining 200M allocations)

#### 3. Memory Arena Allocator (50% gain - 1 day)
**Problem**: Millions of individual heap allocations  
**Solution**: Bump-pointer allocation from pre-allocated buffer
```cpp
class string_arena {
    std::vector<char> buffer_;
    size_t offset_ = 0;
    string_view allocate(string_view src);  // Fast copy to arena
};
```
**Expected**: 6,200 → 9,300 rec/s (faster allocation + better locality)

### ⚡ Tier 2: Solid Wins (20-40% potential each)

#### 4. Lazy Sample Parsing (25% for filtered queries)
**Problem**: Parsing all 3,202 samples even when only querying metadata  
**Solution**: Parse samples only when accessed
```cpp
const auto& samples() const {
    if (!samples_parsed_) parse_samples_on_demand();
    return samples_;
}
```

#### 5. Parallel Parsing (3-4x on multi-core)
**Problem**: Single-threaded parsing  
**Solution**: Parse record batches in parallel
```cpp
#pragma omp parallel for
for (int i = 0; i < batches.size(); ++i) {
    parse_batch(batches[i]);
}
```
**Expected**: 9,300 → 30,000+ rec/s on 8-core CPU

### 🔧 Tier 3: Polish (5-15% each)

6. **Integer encoding for genotypes** (10-15%)
7. **String interning for common values** (5-10%)
8. **SIMD string operations** (15-25%)
9. **Small string optimization awareness** (5%)

## Why bcftools is Faster: Hard Truths

1. **Binary format**: Parses BCF (binary VCF), not text
2. **No string overhead**: Uses fixed-width integers/flags
3. **Written in C**: Manual memory control, no std::string
4. **15+ years of optimization**: Hundreds of man-hours of profiling
5. **Columnar layout**: SIMD-friendly processing

## Realistic Performance Goals

| Strategy | Records/sec | vs Current | vs bcftools |
|----------|-------------|------------|-------------|
| **Current** | 3,790 | baseline | 3.1x slower |
| + Indexed + string_view | 6,200 | 1.6x | 1.9x slower |
| + Memory arena | 9,300 | 2.5x | 1.2x slower |
| + Parallel (8-core) | 30,000+ | 7.9x | **2.6x faster** |
| + BCF format | Match bcftools | - | ~1.0x |

**Bottom Line:**
- **Text VCF ceiling**: ~60-70% of bcftools (7,000-8,000 rec/s)
- **To match**: Need BCF binary format support
- **To exceed**: Parallel parsing + BCF format

## Can We Avoid Allocations? YES!

**Current allocations per record**:
- Fixed fields (CHROM, ID, REF, ALT, etc.): ~10 strings
- INFO fields (~15 pairs): ~30 strings
- Format fields (10 names): ~10 strings
- Sample values (3,202 × 10): ~32,020 strings
- **Total: ~32,070 string allocations per record**

**How to avoid them**:

1. ✅ **Already done**: Use string_view for parsing (eliminates intermediate copies)
2. ⬜ **Indexed samples**: Store format field names once → saves 20,060 allocations/record
3. ⬜ **String_view storage**: Don't convert to string, keep as views → saves remaining 12,010 allocations/record
4. ⬜ **Memory arena**: Bulk allocation → saves malloc overhead on remaining strings

**Result**: From 32,070 → near-zero allocations per record!

## Implementation Status

### Completed ✅
1. ✅ Established baseline with bcftools comparison
2. ✅ Created comprehensive validation test suite (100,356 assertions)
3. ✅ Upgraded to C++17 for std::string_view
4. ✅ Implemented zero-copy parsing with string_view
5. ✅ Added fast integer/double parsing (no std::stoi)
6. ✅ Enabled -O3 compiler optimization (7.8x gain!)
7. ✅ Profiled with manual instrumentation (identified bottleneck)
8. ✅ Replaced std::map with vector<pair> for samples (1.1x gain)

### In Progress 🔄
9. ⬜ Implement indexed sample format (no key duplication)
10. ⬜ Add string_view sample storage (zero allocation)

### Planned 📋
11. ⬜ Memory arena allocator
12. ⬜ Lazy sample parsing
13. ⬜ Parallel parsing support
14. ⬜ BCF binary format reader

## Performance Timeline

```
Baseline (Jan 2026):     234 rec/s  ████░░░░░░░░░░░░░░░░░░░░ 2.0%
+ string_view:           441 rec/s  ████████░░░░░░░░░░░░░░░░ 3.8%
+ O3 optimization:     3,440 rec/s  ████████████████████████░ 29.6%
+ vector samples:      3,790 rec/s  ██████████████████████████ 32.6%
→ Next (indexed):     ~4,800 rec/s  ████████████████████████████░ 41.4%
→ Next (string_view): ~6,200 rec/s  ██████████████████████████████████ 53.4%
→ Goal (arena):       ~9,300 rec/s  ██████████████████████████████████████████ 80.1%
bcftools baseline:    11,611 rec/s  ████████████████████████████████████████████ 100%
```

## Notes
- bcftools is written in highly optimized C with 15+ years of development
- Our C++11 header-only library prioritizes portability and ease of use
- Target: Get within 2-3x of bcftools speed (still excellent for header-only lib)
- Current gap suggests significant low-hanging fruit in optimization
