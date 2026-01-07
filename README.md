# Boost.Genetics

A lightweight, header-only C++ library for working with VCF (Variant Call Format) files, designed as a prototype candidate for the Boost C++ Libraries.

## Features

- **Header-only**: No compilation required, just include and use
- **VCFv4.3 compliant**: Full support for the VCF specification
- **Modern C++**: C++11 standard
- **Well-tested**: Comprehensive test suite using Catch2
- **Error context**: Detailed error messages with line numbers

## VCF Support

Complete implementation of the Variant Call Format specification (VCFv4.3):
- **VCF file reading and writing**
- **Full header metadata support**: INFO, FORMAT, FILTER, contig definitions
- **Sample genotype data handling**: Multiple samples with arbitrary FORMAT fields
- **Standard-compliant parsing**: SNPs, insertions, deletions, and complex variants
- **Error handling**: Custom exceptions with line number context
- **Missing value support**: Proper handling of `.` values

## Quick Start

### Using the Library

Since this is a header-only library, simply include the headers you need:

```cpp
#include <genetics/vcf.hpp>
#include <iostream>

int main() {
    using namespace boost::genetics;
    
    // Read VCF file
    try {
        vcf::reader reader("variants.vcf");
        std::cout << "VCF version: " << reader.version() << std::endl;
        
        vcf::record rec;
        while (reader.read_record(rec)) {
            std::cout << rec.chrom() << ":" << rec.pos() 
                      << " " << rec.ref() << ">" << rec.alt()[0] << std::endl;
        }
    } catch (const vcf::vcf_parse_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "Line: " << e.line() << std::endl;
    }
    
    return 0;
}
```

### Writing VCF Files

```cpp
#include <genetics/vcf.hpp>

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

### Building with CMake

```bash
mkdir build
cd build
cmake ..
make
```

### Running Examples

```bash
# Build examples
cmake -DBUILD_EXAMPLES=ON ..
make

# Run VCF example
./vcf_example
```

### Running Tests

```bash
# Build tests
cmake -DBUILD_TESTS=ON ..
make

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

- ✅ **VCFv4.3 text format** - Fully implemented
- 🔄 **Phase 2 validation** - Type-safe INFO/FORMAT parsing
- 📋 **Structural variants** - Support for `<DEL>`, `<INS>`, breakends
- 📋 **gVCF support** - Reference blocks with `<*>` allele
- 📋 **BCF binary format** - Binary compressed VCF
- 📋 **BGZF compression** - Block compression support
- 📋 **Tabix indexing** - Random access by genomic region
