# VCF Parser Optimization Strategy

## Executive Summary

**Progress Update (Jan 2026)**:
- Baseline: **234 records/sec** 
- Current: **3,790 records/sec** (16.2x improvement ✅)
- bcftools: **11,611 records/sec**
- Gap: **3.1x** (down from 50x!)

**Deep Audit Findings**: Identified **401 million string allocations** as primary bottleneck.  
**Next Steps**: Indexed format + string_view storage → projected 6,200-9,300 rec/s (1.2-1.9x of bcftools)

This document outlines systematic optimization strategies to close the remaining performance gap.

---

## 0. Deep Performance Audit Summary

### Bottleneck Analysis (Profiling Results)

**Time breakdown**:
```
I/O:          0.0002s  (0.01%)  ← Not the bottleneck!
Parsing:      1.6540s  (99.9%)  ← THE BOTTLENECK
Other:        0.0007s  (0.04%)
```

**Within parse_record() (99.9% of time)**:
- Sample parsing: 95% (3,202 samples × 10 fields each)
- INFO parsing: 3%
- Fixed fields: 2%

### Root Cause: String Allocation Tsunami

**Per-record allocations** (6,272 records in test file):
- Fixed fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER): ~10 strings
- INFO fields (~15 key-value pairs): ~30 strings  
- FORMAT field names (10 fields): ~10 strings
- Sample values (3,202 samples × 10 fields): **32,020 strings**
- **Total per record: ~32,070 string allocations**

**Full file parse**: 32,070 × 6,272 = **201,206,880 string objects created**  
**Including keys**: Each sample stores keys too = **401,658,880 total allocations**

### Micro-Benchmark: Storage Strategy Comparison

Testing 3,202 samples × 6,272 records × 10 fields = 200,829,440 values:

| Strategy | Records/sec | Allocations | Speedup |
|----------|-------------|-------------|---------|
| **Current** (pair<string,string>) | 4,293 | 401M | baseline |
| **Indexed** (vector<string>, no keys) | 5,468 | 200M | ⬆ **27%** |
| **String_view** (zero-copy) | 5,575 | 0 | ⬆ **30%** |

**Key Insight**: We're duplicating format field names 3,202 times per record!  
**Solution**: Store format fields once, samples are just value arrays.

---

## 1. Profiling & Measurement

### 1.1 Profiling Tools

**macOS (Primary Platform)**
```bash
# Instruments (Xcode) - Time Profiler
instruments -t "Time Profiler" ./build/benchmark_vcf_parsing
# Export trace, analyze hot paths

# Sample-based profiling
sample ./build/benchmark_vcf_parsing 10 -f profiling_output.txt
```

**Linux (If Available)**
```bash
# perf with call graphs
perf record -g ./build/benchmark_vcf_parsing
perf report --stdio

# Cachegrind for cache misses
valgrind --tool=cachegrind ./build/benchmark_vcf_parsing
cg_annotate cachegrind.out.<pid>
```

**Manual Instrumentation**
```cpp
// Add timing points in vcf.hpp
auto start = std::chrono::high_resolution_clock::now();
// ... parsing code ...
auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
    std::chrono::high_resolution_clock::now() - start);
```

### 1.2 Metrics to Track

- **Records/second** (primary metric)
- **CPU cycles per record**
- **Memory allocations per record**
- **Cache miss rate** (L1/L2/L3)
- **Branch misprediction rate**
- **Peak memory usage**

---

## 2. High-Impact Optimizations

### ✅ COMPLETED: Phase 1 Optimizations (16.2x achieved)

#### 2.1 String View Parsing ✅ DONE
**Impact**: 1.88x improvement (234 → 441 rec/s)  
**Implementation**: Upgraded to C++17, used std::string_view throughout parse_record()
```cpp
// Before: Many allocations
std::string chrom = line.substr(0, tab1);
std::string pos_str = line.substr(tab1+1, tab2-tab1-1);

// After: Zero-copy parsing
std::string_view chrom(line.data(), tab1);
int pos = detail::parse_uint_fast(line.data() + tab1 + 1);
```

#### 2.2 Fast Integer Parsing ✅ DONE
**Impact**: Included in string_view improvements  
**Implementation**: Custom parse_uint_fast() and parse_double_fast()
```cpp
inline unsigned int parse_uint_fast(const char* str) noexcept {
    unsigned int result = 0;
    while (*str >= '0' && *str <= '9') {
        result = result * 10 + (*str - '0');
        ++str;
    }
    return result;
}
```

#### 2.3 Compiler Optimization ✅ DONE  
**Impact**: 7.8x improvement (441 → 3,440 rec/s) - **BIGGEST WIN!**  
**Implementation**: Added -O3 -g -fno-omit-frame-pointer to CMakeLists.txt
```cmake
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -fno-omit-frame-pointer")
```

#### 2.4 Vector-Based Sample Storage ✅ DONE
**Impact**: 1.1x improvement (3,440 → 3,790 rec/s)  
**Implementation**: Replaced std::map<string,string> with vector<pair<string,string>>
```cpp
// Before: Red-black tree with poor cache locality
std::vector<std::map<std::string, std::string>> samples_;

// After: Sequential memory access
std::vector<std::vector<std::pair<std::string, std::string>>> samples_;
```

---

### 🔥 NEXT: Phase 2 - Eliminate Remaining Allocations

#### 2.5 Indexed Sample Format (Expected: 27% gain, 2 hours)

**Problem:** Format field keys duplicated for every sample
```cpp
// Current: "GT" stored 3,202 times per record!
samples_[0] = {{"GT", "0/1"}, {"DP", "50"}, {"GQ", "99"}};
samples_[1] = {{"GT", "1/1"}, {"DP", "45"}, {"GQ", "95"}};  // "GT", "DP", "GQ" duplicated
// ... 3,200 more times
```

**Solution:** Store format fields once per record
```cpp
class record {
    std::vector<std::string> format_fields_;  // Once: ["GT", "DP", "GQ"]
    std::vector<std::vector<std::string>> sample_values_;  // Just values
    
    // Accessor by field name
    std::string get_sample_value(size_t sample_idx, const std::string& field) const {
        for (size_t i = 0; i < format_fields_.size(); ++i) {
            if (format_fields_[i] == field) return sample_values_[sample_idx][i];
        }
        return "";
    }
};
```

**Expected Impact:** 401M → 200M allocations, 3,790 → 4,800 rec/s

---

#### 2.6 String_view Sample Values (Expected: 30% gain, 4 hours)

**Problem:** Converting every sample value to std::string
```cpp
// Current: Allocates 200M+ strings
sample_vec.emplace_back(format_fields[i], std::string(value));  // Allocation!
```

**Solution:** Keep line buffer alive, use string_view
```cpp
class record {
    std::shared_ptr<std::string> line_buffer_;  // Keep source alive
    std::vector<std::vector<std::string_view>> sample_values_;  // No allocations!
    
public:
    void set_line_buffer(std::shared_ptr<std::string> buffer) {
        line_buffer_ = buffer;
    }
};

// In reader::parse_record()
auto line_buffer = std::make_shared<std::string>(std::move(line));
rec.set_line_buffer(line_buffer);

// Parse directly into string_views
std::vector<std::string_view> values;
detail::split_view(sample_view, ':', [&](std::string_view v) {
    values.push_back(v);  // No allocation, points into line_buffer
});
```

**Pros**: Zero allocation, cache-friendly  
**Cons**: Requires buffer lifetime management, can't modify samples in-place  
**Expected Impact:** 200M → 0 allocations, 4,800 → 6,200 rec/s

---

#### 2.7 Memory Arena Allocator (Expected: 50% gain, 1 day)

**Problem:** Even with string_view for samples, still allocating INFO, fixed fields  
**Solution:** Bump-pointer allocation from pre-allocated buffer
```cpp
class string_arena {
    std::vector<char> buffer_;
    size_t offset_ = 0;
    
public:
    string_arena(size_t initial_size = 1024 * 1024) : buffer_(initial_size) {}
    
    std::string_view allocate(std::string_view src) {
        // Grow if needed
        if (offset_ + src.size() > buffer_.size()) {
            buffer_.resize(buffer_.size() * 2);
        }
        
        // Copy to arena
        char* ptr = &buffer_[offset_];
        std::memcpy(ptr, src.data(), src.size());
        offset_ += src.size();
        
        return {ptr, src.size()};
    }
    
    void reset() { offset_ = 0; }  // Bulk deallocation
};

class reader {
    string_arena arena_;
    
public:
    bool parse_record(record& rec) {
        arena_.reset();  // Reset for this record
        
        // Allocate from arena instead of heap
        rec.set_chrom(arena_.allocate(chrom_view));
        rec.set_id(arena_.allocate(id_view));
        // ...
    }
};
```

**Pros**: Fast allocation (just pointer bump), great locality, bulk deallocation  
**Cons**: Requires per-record or per-batch reset  
**Expected Impact:** 6,200 → 9,300 rec/s

---

### 2.3 Optimized String Splitting (Expected: 3-5x improvement)

**Problem:** Repeated `find()` and `substr()` calls
```cpp
// Current: Many allocations
size_t pos = 0;
while ((pos = line.find('\t')) != std::string::npos) {
    fields.push_back(line.substr(0, pos));  // Allocation!
    line = line.substr(pos + 1);            // Allocation!
}
```

**Solution:** Single-pass field extraction
```cpp
struct field_iterator {
    const char* current;
    const char* end;
    char delimiter;
    
    string_view next() {
        const char* start = current;
        while (current != end && *current != delimiter) ++current;
        string_view result(start, current - start);
        if (current != end) ++current;  // Skip delimiter
        return result;
    }
};
```

**Expected Impact:** 3-5x on field extraction

---

### 2.4 INFO Field Parsing Optimization (Expected: 2-4x improvement)

**Problem:** INFO parsing does multiple passes
```cpp
// Current: Parse, split, allocate map
std::map<std::string, std::string> info;  // Red-black tree, slow inserts
for each key=value pair:
    allocate strings
    insert into map  // O(log n) insertion
```

**Solution:** Optimized INFO parsing
```cpp
// Use flat_map or vector for better cache locality
std::vector<std::pair<string_view, string_view>> info_fields;
info_fields.reserve(estimated_count);  // Avoid reallocation

// Single-pass parse
const char* p = info_str;
while (*p) {
    string_view key = read_until(p, '=');
    string_view val = read_until(p, ';');
    info_fields.emplace_back(key, val);
}

// Convert to strings only when accessed
std::string get_info(const char* key) const {
    for (auto& kv : info_fields) {
        if (kv.first == key) return std::string(kv.second.data, kv.second.size);
    }
    return "";
}
```

**Expected Impact:** 2-4x on INFO-heavy files

---

### 2.5 I/O Optimization (Expected: 1.5-2x improvement)

**Problem:** Line-by-line reading with default buffers

**Solutions:**

**A. Larger Buffer Sizes**
```cpp
// Increase buffer from default 8KB to 1MB
std::ifstream file(filename);
constexpr size_t buffer_size = 1024 * 1024;  // 1MB
std::vector<char> buffer(buffer_size);
file.rdbuf()->pubsetbuf(buffer.data(), buffer_size);
```

**B. Memory-Mapped Files**
```cpp
// POSIX mmap (if available)
#include <sys/mman.h>
#include <fcntl.h>

int fd = open(filename, O_RDONLY);
struct stat sb;
fstat(fd, &sb);
char* data = (char*)mmap(nullptr, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);

// Parse directly from mapped memory
const char* end = data + sb.st_size;
for (const char* line_start = data; line_start < end; ) {
    const char* line_end = (const char*)memchr(line_start, '\n', end - line_start);
    parse_record(string_view(line_start, line_end - line_start));
    line_start = line_end + 1;
}

munmap(data, sb.st_size);
```

**C. Parallel I/O & Parsing**
```cpp
// Read chunks in separate thread, parse in main thread
std::thread io_thread([&]() {
    // Read large chunks into ring buffer
});

// Main thread parses from buffer
```

**Expected Impact:** 1.5-2x improvement

---

### 2.6 Small String Optimization (SSO) Awareness

**Analysis:** Many VCF fields are small (<16 chars)
- CHROM: "chr1", "chr22" (4-5 chars)
- REF/ALT: "A", "T", "G", "C" (1 char)
- ID: "rs12345" (7 chars)

**Strategy:** Use std::string wisely, it has SSO for small strings
```cpp
// Good: These stay on stack (no heap allocation)
std::string ref = "A";      // 1 char
std::string alt = "G";      // 1 char
std::string id = "rs1234";  // 6 chars

// Consider: Custom fixed-size string for nucleotides
struct nucleotide {
    char data[4];
    uint8_t len;
    // Optimized for 1-3 char sequences
};
```

---

### 2.7 Container Optimization

**Problem:** `std::map` for INFO/FORMAT has poor cache locality

**Solutions:**

**A. Flat Map (Sorted Vector)**
```cpp
template<typename K, typename V>
class flat_map {
    std::vector<std::pair<K, V>> data_;
    
    void insert(K key, V val) {
        auto it = std::lower_bound(data_.begin(), data_.end(), key);
        data_.insert(it, {key, val});
    }
    
    V* find(const K& key) {
        auto it = std::lower_bound(data_.begin(), data_.end(), key);
        return (it != data_.end() && it->first == key) ? &it->second : nullptr;
    }
};
```

**B. Small Vector Optimization**
```cpp
// Most records have 0-5 ALT alleles
template<typename T, size_t N>
class small_vector {
    T stack_storage[N];
    T* heap_storage = nullptr;
    size_t size_ = 0;
    // Use stack storage for size <= N
};

small_vector<std::string, 4> alt_alleles;  // No heap for ≤4 alleles
```

**Expected Impact:** 1.5-2x for INFO/sample-heavy files

---

### 2.8 Branch Prediction & Loop Optimization

**Problem:** Unpredictable branches in hot loops

**Solutions:**

**A. Branchless Parsing**
```cpp
// Instead of:
if (c == '\t') { /* ... */ }
else if (c == '\n') { /* ... */ }
else { /* ... */ }

// Use lookup table:
enum char_type { REGULAR, TAB, NEWLINE };
char_type type_table[256] = { /* ... */ };

switch (type_table[static_cast<unsigned char>(c)]) {
    case TAB: /* ... */ break;
    case NEWLINE: /* ... */ break;
    case REGULAR: /* ... */ break;
}
```

**B. Loop Unrolling**
```cpp
// Parse 4 characters at a time (SIMD-style)
while (remaining >= 4) {
    // Process 4 bytes
    remaining -= 4;
}
// Handle remainder
```

**Expected Impact:** 1.2-1.5x improvement

---

### 2.9 Validation Optimization

**Current:** Full validation on every record (expensive)

**Strategy:** Lazy/Optional Validation
```cpp
enum class validation_level {
    none,        // No validation (fastest)
    basic,       // Critical checks only (REF matches [ACGTN]+)
    standard,    // Standard checks
    strict       // Full VCF v4.3 compliance
};

reader.set_validation(validation_level::basic);  // User choice
```

**Expected Impact:** 2-3x when validation disabled

---

## 3. Implementation Roadmap (Updated)

### ✅ Phase 1: Foundation & String Optimization (COMPLETED)
**Goal:** Establish profiling baseline and eliminate obvious inefficiencies

1. ✅ Upgraded to C++17 for std::string_view
2. ✅ Implemented zero-copy parsing with string_view
3. ✅ Added fast integer/double parsers
4. ✅ Enabled -O3 compiler optimization
5. ✅ Added profiling instrumentation
6. ✅ Replaced std::map with vector<pair> for samples

**Result:** 234 → 3,790 rec/s (16.2x improvement)  
**Deliverable:** Validated optimizations with 100,356 assertions vs bcftools ✅

---

### 🔄 Phase 2: Eliminate Allocations (IN PROGRESS - 1 week)
**Goal:** Remove remaining 200M+ string allocations

1. ⬜ **Indexed sample format** (2 hours)
   - Remove format field keys from sample vectors
   - Store format_fields_ once per record
   - Update API to access by index or field name

2. ⬜ **String_view sample storage** (4 hours)
   - Add line_buffer_ to record class
   - Store string_views instead of strings for sample values
   - Test memory lifetime carefully

3. ⬜ **Benchmark and validate** (2 hours)
   - Run full test suite
   - Verify 100,356 assertions still pass
   - Measure performance gain

**Expected:** 3,790 → 6,200 rec/s (1.6x improvement)  
**Deliverable:** Near-zero allocation parser

---

### 📋 Phase 3: Memory Arena (PLANNED - 3 days)
**Goal:** Optimize remaining allocations with arena allocator

1. ⬜ Implement string_arena class
2. ⬜ Integrate with reader for INFO/fixed fields  
3. ⬜ Add per-record or per-batch reset
4. ⬜ Benchmark and tune arena size

**Expected:** 6,200 → 9,300 rec/s (1.5x improvement)  
**Deliverable:** Production-ready high-performance parser

---

### 📋 Phase 4: Advanced Features (FUTURE - 1-2 weeks)

1. ⬜ **Lazy sample parsing** (20-30% for filtered queries)
   - Don't parse samples unless accessed
   - Useful for queries that only need metadata

2. ⬜ **Parallel parsing** (3-4x on multi-core)
   - Parse records in parallel batches
   - OpenMP or std::thread implementation

3. ⬜ **BCF binary format reader**
   - Match bcftools directly
   - Requires BGZF decompression

**Expected:** 9,300 → 30,000+ rec/s (parallel) or match bcftools (BCF)  
**Deliverable:** Feature-complete high-performance library

---

## 4. Expected Performance Trajectory

**Actual Progress:**
```
Phase 1 Complete:
  Baseline:          234 rec/s  ██░░░░░░░░░░░░░░░░░░░░ 2.0%
  + string_view:     441 rec/s  ███░░░░░░░░░░░░░░░░░░░ 3.8%
  + O3:            3,440 rec/s  ████████████████████░░░ 29.6%
  + vector:        3,790 rec/s  ██████████████████████░ 32.6%
```

**Projected Next Steps:**
```
Phase 2: Zero Allocation
  + Indexed:      ~4,800 rec/s  ████████████████████████░ 41.3%
  + String_view:  ~6,200 rec/s  ███████████████████████████░ 53.4%

Phase 3: Memory Arena  
  + Arena:        ~9,300 rec/s  ████████████████████████████████░ 80.1%

Phase 4: Advanced
  + Parallel:    ~30,000 rec/s  ██████████████████████████████████████████ 258%
  
bcftools baseline: 11,611 rec/s  ████████████████████████████████████ 100%
```

**Key Milestones:**
- ✅ **16.2x from baseline** (234 → 3,790 rec/s)
- ⏭️ **Within 2x of bcftools** (~6,200 rec/s with indexed + string_view)
- 🎯 **Match bcftools** (~9,300 rec/s with arena, or BCF support)
- 🚀 **Exceed bcftools** (~30,000 rec/s with parallel parsing)

---

## 5. Optimization Guidelines

### 5.1 Principles

1. **Measure First:** Profile before optimizing
2. **One Change at a Time:** Isolate impact
3. **Preserve Correctness:** All tests must pass
4. **Maintain Readability:** Header-only library needs clarity
5. **Document Trade-offs:** Note when accuracy/safety traded for speed

### 5.2 Benchmarking Process

```bash
# Before each optimization
./build/benchmark_vcf_parsing > baseline.txt

# After optimization
./build/benchmark_vcf_parsing > optimized.txt

# Compare
diff baseline.txt optimized.txt

# Run validation tests
./build/test_validation

# Check no regressions
diff tests/data/reference/bcftools_output.txt tests/data/reference/our_output.txt
```

### 5.3 Safety Checklist

- [ ] All unit tests pass (41 tests, 381 assertions)
- [ ] All validation tests pass (100,356 assertions vs bcftools)
- [ ] No memory leaks (valgrind clean)
- [ ] No undefined behavior (sanitizers clean)
- [ ] Performance improvement documented
- [ ] Code review completed

---

## 6. Advanced Techniques (Future)

### 6.1 SIMD Optimization

- Use SSE/AVX for string operations
- Vectorized integer parsing
- Parallel validation checks

**Expected:** Additional 1.5-2x on modern CPUs

### 6.2 Compile-Time Optimization

```cpp
// Template specialization for common cases
template<> class record_parser<validation_level::none> {
    // Fastest path, no checks
};

template<> class record_parser<validation_level::strict> {
    // Full validation
};
```

### 6.3 Just-In-Time Compilation

- Compile parser for specific VCF header
- Optimize field extraction for known formats
- Cache compiled parsers

**Expected:** 2-3x for repeated parsing of similar files

---

## 7. Comparison with bcftools

### Why bcftools is faster:

1. **C implementation:** Lower-level, manual memory management
2. **BCF focus:** Binary format optimized for speed
3. **Decades of optimization:** HTSlib highly tuned
4. **SIMD usage:** Uses vectorization
5. **Custom allocators:** Pool allocators for common objects

### Our Advantages:

1. **Header-only:** No linking, easy integration
2. **C++17:** Modern API, std::string_view, if constexpr, structured bindings
3. **Readable code:** Easier to maintain/extend
4. **Template flexibility:** Generic programming benefits

### Realistic Target:

**Within 2-3x is excellent** for:
- Header-only library
- C++17 standard (native string_view, constexpr improvements)
- Emphasis on code clarity
- Full validation support

---

## 8. Monitoring & Metrics

### Performance Dashboard (Updated Jan 2026)

```
═══════════════════════════════════════════════════════════
             VCF PARSER PERFORMANCE STATUS
═══════════════════════════════════════════════════════════

Parse Speed:      3,790 rec/sec  ██████████████████████░░ 32.6%
Target (2x):      5,805 rec/sec  ███████████████████████░ 50.0%
Target (1.2x):    9,675 rec/sec  ████████████████████████ 83.3%
bcftools:        11,611 rec/sec  ████████████████████████ 100%

Progress:         16.2x from baseline ✅
Gap to bcftools:  3.1x (down from 50x!)
Phase:            Phase 1 Complete, Phase 2 Next

Allocations:      ~32,000 per record
Next Goal:        Reduce to near-zero

Memory Usage:     ~13.9 KB/record
Validation:       ✅ 100% match (100,356 assertions)
Test Coverage:    ✅ All tests passing

Recent Optimizations:
  ✅ C++17 upgrade (string_view support)
  ✅ Zero-copy field parsing
  ✅ Fast integer/double parsers
  ✅ -O3 compiler optimization (7.8x gain!)
  ✅ Vector-based sample storage (1.1x gain)

Next Up:
  ⏭️ Indexed sample format (27% projected)
  ⏭️ String_view sample storage (30% projected)
═══════════════════════════════════════════════════════════
```

### Allocation Tracking

```
Current allocations per record: ~32,070
  - Fixed fields:     ~10
  - INFO fields:      ~30
  - Format fields:    ~10
  - Sample values:    ~32,020  ← 99.8% of allocations!

After indexed format: ~12,010 (62.5% reduction)
After string_view:    ~10-20 (99.9% reduction)
```

---

## References

- **VCF Spec:** https://samtools.github.io/hts-specs/VCFv4.3.pdf
- **bcftools source:** https://github.com/samtools/bcftools
- **HTSlib source:** https://github.com/samtools/htslib
- **C++ optimization guides:** 
  - Agner Fog's optimization manuals
  - "Optimized C++" by Kurt Guntheroth
  - CppCon performance talks

---

## Conclusion

With systematic optimization focusing on:
1. Zero-copy string parsing (std::string_view)
2. Fast integer conversion
3. Optimized containers
4. Efficient I/O

We can realistically achieve **7,000-9,500 records/sec**, placing us within 1.2-1.7x of bcftools while maintaining:
- Clean, maintainable code
- Full VCF v4.3 compliance
- Header-only convenience
- C++17 modern features

**Next Step:** Begin Phase 1 profiling to validate assumptions and identify actual hotspots.
