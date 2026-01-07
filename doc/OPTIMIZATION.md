# VCF Parser Optimization Strategy

## Executive Summary

Current baseline: **234 records/sec** vs bcftools: **11,826 records/sec** (50x gap)  
Target: **3,900-5,900 records/sec** (within 2-3x of bcftools - excellent for header-only C++11 library)

This document outlines systematic optimization strategies to close the performance gap.

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

### 2.1 String View Parsing (Expected: 10-30x improvement)

**Problem:** Every field creates new `std::string` allocations
```cpp
// Current: Allocates 10+ strings per record
std::string chrom = line.substr(0, tab1);
std::string pos_str = line.substr(tab1+1, tab2-tab1-1);
int pos = std::stoi(pos_str);  // Another allocation in stoi
```

**Solution:** Zero-copy parsing with string_view
```cpp
// C++17: Use std::string_view natively
#include <string_view>

// Parse without allocation
std::string_view chrom(line.data(), tab1);
int pos = parse_int_fast(line.data() + tab1 + 1, tab2 - tab1 - 1);
```

**Implementation Plan:**
1. Use `std::string_view` from C++17 standard library
2. Add fast integer parsers (avoid `std::stoi`, `std::strtol`)
3. Rewrite `parse_line()` to use views
4. Convert to `std::string` only when storing

**Expected Impact:** 15-25x speedup (most allocations eliminated)

---

### 2.2 Fast Integer Parsing (Expected: 2-3x improvement)

**Problem:** `std::stoi()` is slow (handles locales, exceptions, edge cases)

**Solution:** Custom integer parser
```cpp
inline int parse_int_fast(const char* str, size_t len) noexcept {
    int result = 0;
    size_t i = 0;
    bool negative = false;
    
    if (len > 0 && str[0] == '-') {
        negative = true;
        i = 1;
    }
    
    for (; i < len && str[i] >= '0' && str[i] <= '9'; ++i) {
        result = result * 10 + (str[i] - '0');
    }
    
    return negative ? -result : result;
}

// For quality scores (float)
inline double parse_double_fast(const char* str, size_t len) noexcept;
```

**Expected Impact:** 2-3x on position/quality parsing

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

## 3. Implementation Roadmap

### Phase 1: Foundation (Week 1)
**Goal:** Establish profiling baseline

1. ✅ Run baseline benchmarks
2. Profile with Instruments/perf
3. Identify top 3 hotspots
4. Document findings

**Deliverable:** Profiling report with hotspot analysis

---

### Phase 2: String Optimization (Week 2)
**Goal:** Eliminate string allocations

1. Implement `string_view` wrapper
2. Rewrite field parsing to use views
3. Add fast integer/float parsers
4. Benchmark improvements

**Expected:** 15-20x improvement  
**Deliverable:** Optimized parsing core

---

### Phase 3: Container Optimization (Week 3)
**Goal:** Improve data structure performance

1. Implement flat_map for INFO/FORMAT
2. Add small_vector for ALT alleles
3. Optimize sample data storage
4. Benchmark improvements

**Expected:** Additional 1.5-2x improvement  
**Deliverable:** Optimized containers

---

### Phase 4: I/O Optimization (Week 4)
**Goal:** Maximize I/O throughput

1. Add large buffer option
2. Implement memory-mapped file support
3. Consider parallel parsing
4. Benchmark improvements

**Expected:** Additional 1.5-2x improvement  
**Deliverable:** Optimized I/O layer

---

### Phase 5: Polish & Tuning (Week 5)
**Goal:** Final optimizations

1. Profile again, find remaining hotspots
2. Implement micro-optimizations
3. Add optional validation levels
4. Final benchmarks

**Expected:** Additional 1.2-1.5x improvement  
**Deliverable:** Production-ready optimized parser

---

## 4. Expected Final Performance

**Conservative Estimate:**
- Phase 2 (strings): 15x → **3,510 records/sec**
- Phase 3 (containers): 1.5x → **5,265 records/sec**
- Phase 4 (I/O): 1.5x → **7,898 records/sec**
- Phase 5 (polish): 1.2x → **9,478 records/sec**

**Final:** ~9,500 records/sec (**1.25x slower than bcftools**)  
✅ **Exceeds 2-3x target!**

**Optimistic Estimate:**
- Phases combined could reach 25-30x improvement
- **~7,000-9,000 records/sec range likely**

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

### Performance Dashboard

```
Current Status: Phase 1 (Baseline)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Parse Speed:        234 rec/sec  [████░░░░░░░░░░░░░░░░] 2.0%
Target (2x):      5,913 rec/sec
Target (3x):      3,942 rec/sec
bcftools:        11,826 rec/sec

Memory Usage:     ~1.2 MB/1000 records
Validation:       ✅ 100% match with bcftools
Test Coverage:    ✅ 100,356 assertions passing
```

### Update After Each Phase

Track progress toward target with each optimization phase.

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
