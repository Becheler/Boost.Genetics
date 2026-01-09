#!/bin/bash
# Comprehensive benchmark runner for Boost.Genetics
# Runs memory, scalability, and writer benchmarks

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if data file exists
DATA_FILE="${1:-tests/data/reference/chr22_subset.vcf}"
if [ ! -f "$DATA_FILE" ]; then
    echo -e "${RED}Error: Data file not found: $DATA_FILE${NC}"
    echo "Usage: $0 [vcf_file]"
    exit 1
fi

echo -e "${GREEN}=== Boost.Genetics Comprehensive Benchmark Suite ===${NC}"
echo "Data file: $DATA_FILE"
echo ""

# Create output directory
RESULTS_DIR="benchmark_results"
mkdir -p "$RESULTS_DIR"

# Check if binaries exist
if [ ! -f "build/benchmark_memory" ]; then
    echo -e "${RED}Error: Benchmarks not built. Run: cmake -B build -DBUILD_BENCHMARKS=ON && cmake --build build${NC}"
    exit 1
fi

# 1. Memory Profiling Benchmark
echo -e "${YELLOW}[1/4] Running Memory Profiling Benchmark...${NC}"
./build/benchmark_memory "$DATA_FILE" 4 | tee "$RESULTS_DIR/memory_profile.txt"
echo ""

# 2. Writer Performance Benchmark
echo -e "${YELLOW}[2/4] Running Writer Performance Benchmark...${NC}"
./build/benchmark_writer "$DATA_FILE" "$RESULTS_DIR/output_test.vcf" | tee "$RESULTS_DIR/writer_performance.txt"
rm -f "$RESULTS_DIR/output_test.vcf"
echo ""

# 3. Scalability Benchmark (quick mode)
echo -e "${YELLOW}[3/4] Running Scalability Benchmark (quick mode)...${NC}"
./build/benchmark_scalability --quick | tee "$RESULTS_DIR/scalability.txt"
echo ""

# 4. Parallel Scalability Test
echo -e "${YELLOW}[4/4] Running Parallel Scalability Test...${NC}"
./build/benchmark_scalability --parallel "$DATA_FILE" | tee "$RESULTS_DIR/parallel_scaling.txt"
echo ""

# Generate summary report
echo -e "${GREEN}=== Benchmark Summary ===${NC}"
echo ""
echo "Memory Profile:"
grep -A 3 "Peak RSS:" "$RESULTS_DIR/memory_profile.txt" | tail -4 || echo "No memory data"
echo ""
echo "Writer Performance:"
grep "Write rate:" "$RESULTS_DIR/writer_performance.txt" | head -1 || echo "No writer data"
echo ""
echo "Scalability:"
grep "✓" "$RESULTS_DIR/scalability.txt" || grep "⚠" "$RESULTS_DIR/scalability.txt" || echo "No scalability data"
echo ""

echo -e "${GREEN}Results saved to: $RESULTS_DIR/${NC}"
echo "Files:"
ls -lh "$RESULTS_DIR"/*.txt

# Optionally generate markdown report
if command -v pandoc &> /dev/null; then
    echo ""
    echo "Generating markdown report..."
    cat > "$RESULTS_DIR/BENCHMARK_REPORT.md" << EOF
# Boost.Genetics Benchmark Report
Generated: $(date)

## Test Configuration
- Data file: \`$DATA_FILE\`
- System: $(uname -s) $(uname -r)
- CPU: $(sysctl -n machdep.cpu.brand_string 2>/dev/null || grep "model name" /proc/cpuinfo | head -1 | cut -d: -f2)

## Memory Profiling
\`\`\`
$(cat "$RESULTS_DIR/memory_profile.txt")
\`\`\`

## Writer Performance
\`\`\`
$(cat "$RESULTS_DIR/writer_performance.txt")
\`\`\`

## Scalability Analysis
\`\`\`
$(cat "$RESULTS_DIR/scalability.txt")
\`\`\`

## Parallel Scaling
\`\`\`
$(cat "$RESULTS_DIR/parallel_scaling.txt")
\`\`\`
EOF
    echo -e "${GREEN}Markdown report: $RESULTS_DIR/BENCHMARK_REPORT.md${NC}"
fi
