# Boost.Genetics

A lightweight, header-only C++ library for genetics datasets, inspired by Boost.Bloom. This library is designed as a prototype candidate for the Boost C++ Libraries, providing modular and efficient tools for working with genetic data.

## Features

- **Header-only**: No compilation required, just include and use
- **Lightweight**: Minimal dependencies, no heavy tooling required
- **Modular design**: Clear separation of concerns
- **Modern C++**: C++11 standard
- **Well-tested**: Comprehensive test suite using Catch2
- **Benchmarked**: Performance benchmarks using Google Benchmark

## Components

### Sequence (`genetics/sequence.hpp`)
Represents genetic sequences (DNA, RNA, or protein) with utilities for:
- GC content calculation
- Reverse complement generation
- Sequence validation

### Variant (`genetics/variant.hpp`)
Represents genetic variants (SNPs, insertions, deletions) with utilities for:
- Variant type detection
- Transition/transversion classification
- Chromosome and position tracking

## Quick Start

### Using the Library

Since this is a header-only library, simply include the headers you need:

```cpp
#include <genetics/genetics.hpp>
#include <iostream>

int main() {
    using namespace boost::genetics;
    
    // Create a DNA sequence
    sequence dna("ATCGATCG", sequence::type::DNA);
    std::cout << "GC Content: " << dna.gc_content() * 100 << "%" << std::endl;
    
    // Get reverse complement
    sequence rc = dna.reverse_complement();
    std::cout << "Reverse Complement: " << rc.data() << std::endl;
    
    // Create a variant
    variant snp("chr1", 12345, "A", "G");
    std::cout << "Is transition: " << snp.is_transition() << std::endl;
    
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

# Run examples
./sequence_example
./variant_example
```

### Running Tests

```bash
# Build tests
cmake -DBUILD_TESTS=ON ..
make

# Run tests
ctest --output-on-failure
# Or run individual tests
./test_sequence
./test_variant
```

### Running Benchmarks

```bash
# Build benchmarks
cmake -DBUILD_BENCHMARKS=ON ..
make

# Run benchmarks
./benchmark_sequence
./benchmark_variant
```

## Directory Structure

```
Boost.Genetics/
├── include/genetics/     # Header-only library files
│   ├── genetics.hpp      # Main header
│   ├── version.hpp       # Version information
│   ├── sequence.hpp      # Sequence class
│   └── variant.hpp       # Variant class
├── examples/             # Example programs
│   ├── sequence_example.cpp
│   └── variant_example.cpp
├── tests/                # Unit tests (Catch2)
│   ├── catch2/           # Catch2 header
│   ├── test_sequence.cpp
│   └── test_variant.cpp
├── benchmarks/           # Performance benchmarks (Google Benchmark)
│   ├── benchmark_sequence.cpp
│   └── benchmark_variant.cpp
├── cmake/                # CMake configuration files
└── CMakeLists.txt        # Build configuration
```

## Requirements

- C++11 or later
- CMake 3.14 or later (for building examples, tests, and benchmarks)
- Catch2 (included as single-header in `tests/catch2/`)
- Google Benchmark (optional, for benchmarks, will be fetched automatically)

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

Future features may include:
- FASTA/FASTQ file readers
- VCF file parsers
- Sequence alignment utilities
- Population genetics statistics
- More comprehensive variant annotations
