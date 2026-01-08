# Boost.Genetics

[![CI](https://github.com/arnaudbecheler/Boost.Genetics/actions/workflows/ci.yml/badge.svg)](https://github.com/arnaudbecheler/Boost.Genetics/actions/workflows/ci.yml)
[![Benchmark](https://github.com/arnaudbecheler/Boost.Genetics/actions/workflows/benchmark.yml/badge.svg)](https://github.com/arnaudbecheler/Boost.Genetics/actions/workflows/benchmark.yml)
[![Documentation](https://github.com/arnaudbecheler/Boost.Genetics/actions/workflows/pages.yml/badge.svg)](https://arnaudbecheler.github.io/Boost.Genetics/)

A high-performance, header-only C++17 library for genomic data processing, designed as a prototype candidate for the Boost C++ Libraries.

## Performance

**60% faster than bcftools** with parallel parsing:

| Implementation | Records/sec | Speedup |
|----------------|-------------|---------|
| **Boost.Genetics Parallel (8 threads)** | **~19,300** | **~4.5x sequential** |
| Boost.Genetics Sequential | ~4,300 | 1.0x |
| bcftools | ~12,000 | baseline |

*Benchmarks on 1000 Genomes chr22 dataset. See [benchmarking guide](docs/BENCHMARKING.md) for details.*

## Features

- ⚡ **Ultra-fast VCF parsing**: 60% faster than bcftools with parallel processing
- 🔄 **Parallel processing**: Near-linear scaling up to 4 cores
- 💾 **Zero-copy design**: Minimal memory allocations using `string_view`
- 🎯 **Modern C++17**: Clean, type-safe API
- 📦 **Header-only**: No compilation required, just include and use
- ✅ **VCFv4.3 compliant**: Full support for the VCF specification
- 🧪 **Well-tested**: Comprehensive test suite (381 assertions)

## VCF Support

Complete implementation of the Variant Call Format specification (VCFv4.3):

### Phase 1 (Complete)
- **VCF file reading and writing**
- **Full header metadata support**: INFO, FORMAT, FILTER, contig definitions
- **Sample genotype data handling**: Multiple samples with arbitrary FORMAT fields
- **Standard-compliant parsing**: SNPs, insertions, deletions, and complex variants
- **Error handling**: Custom exceptions with line number context
- **Missing value support**: Proper handling of `.` values

### Phase 2 (Complete) 
- **Structured metadata parsing**: Parse ##INFO, ##FORMAT, ##FILTER, ##contig, ##ALT headers into typed structures
- **Percent encoding/decoding**: Full support for encoding special characters in INFO/FORMAT fields
- **Record validation**: Validate INFO fields against header metadata
- **Variant type detection**: Methods to identify SNPs, insertions, deletions
- **Comprehensive metadata access**: Retrieve parsed header definitions programmatically

### Phase 3 (Complete)
- **Structural variant support**: Full support for symbolic alleles (`<DEL>`, `<INS>`, `<DUP>`, `<INV>`, `<CNV>`)
- **Breakend notation**: Complete parser for all 4 breakend patterns (t[p[, t]p], ]p]t, [p[t)
- **SV INFO fields**: Parsing of SVTYPE, END, SVLEN, CIPOS, CIEND, IMPRECISE, MATEID, EVENT
- **SV detection methods**: is_symbolic_sv(), is_breakend(), get_sv_type(), get_sv_info()
- **Confidence intervals**: Support for imprecise SVs with position uncertainty

### Phase 4 (Complete)
- **gVCF reference blocks**: Full support for `<*>` and `<NON_REF>` symbolic alleles
- **END field parsing**: Parse END INFO field for reference block ranges
- **MIN_DP support**: Minimum depth tracking across reference blocks
- **Reference block validation**: Validate genotypes, depth consistency, position ranges
- **Block length calculation**: Automatic calculation of reference block spans
- **Mixed variants**: Support for variants with both alternate alleles and `<*>`

### Validation (Complete)
- **ID pattern validation**: Enforce VCF v4.3 ID naming rules (alphanumeric, underscore, period)
- **Header validation**: Check fileformat line, unique IDs, unique sample names
- **Field validation**: Validate REF (nucleotides), ALT (nucleotides/symbolic/breakend), POS (integers)
- **QUAL validation**: Ensure QUAL is number or dot
- **Metadata validation**: Validate INFO/FORMAT/FILTER keys against header definitions
- **Record validation**: Comprehensive record validation with detailed error reporting
- **Validation result API**: Accumulate and query validation errors

### BCF Binary Format (Complete)
- **BCF v2.2 header parsing**: Read/write BCF magic bytes, version, VCF header text
- **Type encoding/decoding**: Support for INT8, INT16, INT32, FLOAT, CHAR types
- **Type descriptors**: Size and type encoding with overflow handling
- **Genotype encoding**: VCF genotype string (e.g., "0/1", "1|0") to BCF binary format
- **Special values**: MISSING, END_OF_VECTOR markers for integers and floats
- **BGZF compression**: Optional Boost.IOStreams integration for compressed BCF files
- **Conditional compilation**: Automatically enables BCF file I/O if Boost.IOStreams + zlib available

## Quick Start

### Sequential Parsing

```cpp
#include <boost/genetics/vcf.hpp>
#include <iostream>

int main() {
    using namespace boost::genetics;
    
    // Read VCF file
    vcf::reader reader("variants.vcf");
    std::cout << "VCF version: " << reader.version() << std::endl;
    
    vcf::record rec;
    while (reader.read_record(rec)) {
        std::cout << rec.chrom() << ":" << rec.pos() 
                  << " " << rec.ref() << ">" << rec.alt()[0] << std::endl;
    }
    
    return 0;
}
```

### Parallel Parsing (4.5x faster)

```cpp
#include <boost/genetics/vcf.hpp>
#include <iostream>

int main() {
    using namespace boost::genetics;
    
    // Parallel parsing with 8 threads
    vcf::parallel_reader reader("variants.vcf", 8);
    
    // Read header
    auto header = reader.read_header();
    std::cout << "VCF version: " << header.version() << std::endl;
    
    // Parse all records in parallel (unordered for max performance)
    auto records = reader.read_all_unordered();
    
    std::cout << "Parsed " << records.size() << " records" << std::endl;
    for (const auto& rec : records) {
        std::cout << rec.chrom() << ":" << rec.pos() << std::endl;
    }
    
    return 0;
}
```

### Writing VCF Files

```cpp
#include <boost/genetics/vcf.hpp>

int main() {
    using namespace boost::genetics;
    
    // Create writer
    vcf::writer writer("output.vcf", "VCFv4.3");
    
    // Add metadata
    writer.add_contig("chr1", 248956422);
    writer.add_info("DP", "1", "Integer", "Total Depth");
    writer.add_format("GT", "1", "String", "Genotype");
    
    std::vector<std::string> samples = {"Sample1", "Sample2"};
    writer.set_sample_names(samples);
    writer.write_header();
    
    // Create and write variant
    std::vector<std::string> alt = {"G"};
    std::map<std::string, std::string> info;
    info["DP"] = "100";
    
    vcf::record rec("chr1", 12345, "rs123", "A", alt, 30.0, "PASS", info);
    rec.set_format("GT");
    
    std::map<std::string, std::string> sample_data;
    sample_data["GT"] = "0/1";
    rec.add_sample(sample_data);
    rec.add_sample(sample_data);
    
    writer.write_record(rec);
    
    return 0;
}
```

## Installation

### Header-Only (Recommended)

```bash
# Clone the repository
git clone https://github.com/arnaudbecheler/Boost.Genetics.git
cd Boost.Genetics

# Add include directory to your project
# In your CMakeLists.txt:
target_include_directories(your_target PRIVATE ${CMAKE_SOURCE_DIR}/Boost.Genetics/include)
```

### Building with CMake

```bash
mkdir build
cd build
cmake .. -DBUILD_EXAMPLES=ON -DBUILD_TESTS=ON -DBUILD_BENCHMARKS=ON
make -j$(nproc)
```

### Running Examples

```bash
# Sequential parsing
./build/vcf_example

# Parallel parsing benchmark
./build/benchmark_vcf_parallel tests/data/reference/chr22_subset.vcf 8
```

### Running Tests

```bash
# Build and run tests
cmake -DBUILD_TESTS=ON ..
make -j$(nproc)
./build/test_vcf

# Run tests
ctest --output-on-failure
# Or run directly
./test_vcf
```

## Directory Structure

```
Boost.Genetics/
├── include/genetics/     # Header-only library files
│   ├── genetics.hpp      # Main header (includes vcf.hpp)
│   ├── version.hpp       # Version information
│   └── vcf.hpp           # VCF format support
├── examples/             # Example programs
│   └── vcf_example.cpp
├── tests/                # Unit tests (Catch2)
│   ├── catch2/           # Catch2 header
│   └── test_vcf.cpp
├── doc/                  # Documentation
│   └── assets/           # Specification documents
├── cmake/                # CMake configuration files
└── CMakeLists.txt        # Build configuration
```

## Requirements

- C++11 or later
- CMake 3.14 or later (for building examples and tests)
- Catch2 (included as single-header in `tests/catch2/`)

### Optional Dependencies (for BCF support)
- Boost.IOStreams (for BGZF compression)
- zlib (for gzip compression)

If Boost.IOStreams and zlib are found, BCF file I/O is automatically enabled. Otherwise, the library provides BCF encoding/decoding utilities but not file I/O.

## Installation

As a header-only library, you can simply copy the `include/genetics/` directory to your project or install system-wide:

```bash
mkdir build && cd build
cmake ..
sudo make install
```

## License

Distributed under the Boost Software License, Version 1.0.
See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt

## Contributing

This is a prototype candidate library for Boost. Contributions are welcome!

## Roadmap

- ✅ **Phase 1: VCFv4.3 text format** - Core reader/writer
- ✅ **Phase 2: Extended features** - Metadata parsing, validation, percent encoding
- ✅ **Phase 3: Structural variants** - Support for `<DEL>`, `<INS>`, breakends, SVTYPE/END/SVLEN
- ✅ **Phase 4: gVCF support** - Reference blocks with `<*>` allele, END field, MIN_DP
- ✅ **Validation** - Comprehensive VCF v4.3 spec compliance validation
- ✅ **BCF binary format** - Type encoding, genotype encoding, BGZF compression (optional)
- 📋 **Performance optimizations** - Memory pooling, lazy parsing, parallel processing
- 📋 **Advanced features** - Tabix indexing, region queries, VCF merging
