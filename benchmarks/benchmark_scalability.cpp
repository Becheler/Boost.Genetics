//  Scalability benchmark for large VCF files
//  Tests performance with 10K, 100K, 1M, 10M record datasets

#include <boost/genetics/vcf.hpp>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>
#include <cmath>

using namespace boost::genetics;

// Generate a synthetic VCF file with specified number of records
void generate_vcf_file(const std::string& filename, size_t num_records) {
    std::ofstream out(filename);
    
    // Write header
    out << "##fileformat=VCFv4.3\n";
    out << "##contig=<ID=chr1,length=248956422>\n";
    out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
    out << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n";
    out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    out << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n";
    
    // Write records
    for (size_t i = 0; i < num_records; ++i) {
        out << "chr1\t" << (10000 + i) << "\t.\tA\tG\t30\tPASS\t";
        out << "DP=" << (50 + (i % 100)) << ";AF=0." << (i % 10) << "\t";
        out << "GT:DP\t0/1:" << (20 + (i % 50)) << "\t1/1:" << (30 + (i % 50)) << "\n";
    }
}

// Benchmark parsing at different scales
struct benchmark_result {
    size_t num_records;
    double time_seconds;
    double records_per_sec;
    size_t memory_mb;
};

benchmark_result benchmark_sequential(const std::string& filename, size_t expected_records) {
    auto start = std::chrono::high_resolution_clock::now();
    
    vcf::reader reader(filename);
    vcf::record rec;
    size_t count = 0;
    
    while (reader.read_record(rec)) {
        count++;
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    benchmark_result result;
    result.num_records = count;
    result.time_seconds = elapsed.count();
    result.records_per_sec = count / elapsed.count();
    result.memory_mb = 0; // Could add RSS measurement
    
    return result;
}

benchmark_result benchmark_parallel(const std::string& filename, size_t num_threads) {
    auto start = std::chrono::high_resolution_clock::now();
    
    vcf::parallel_reader reader(filename, num_threads);
    auto records = reader.read_all_unordered();
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    benchmark_result result;
    result.num_records = records.size();
    result.time_seconds = elapsed.count();
    result.records_per_sec = records.size() / elapsed.count();
    result.memory_mb = 0;
    
    return result;
}

void print_result(const std::string& name, const benchmark_result& result) {
    std::cout << name << ":\n";
    std::cout << "  Records: " << result.num_records << "\n";
    std::cout << "  Time: " << result.time_seconds << " sec\n";
    std::cout << "  Rate: " << static_cast<int>(result.records_per_sec) << " rec/sec\n";
}

// Test scalability by measuring performance at different dataset sizes
void test_scalability(const std::vector<size_t>& dataset_sizes) {
    std::cout << "\n=== Scalability Test ===" << std::endl;
    std::cout << "Testing with dataset sizes: ";
    for (size_t size : dataset_sizes) {
        std::cout << size << " ";
    }
    std::cout << "\n" << std::endl;
    
    std::vector<benchmark_result> seq_results;
    std::vector<benchmark_result> par_results;
    
    for (size_t size : dataset_sizes) {
        std::string filename = "test_" + std::to_string(size) + ".vcf";
        
        std::cout << "Generating " << size << " records..." << std::flush;
        generate_vcf_file(filename, size);
        std::cout << " done" << std::endl;
        
        std::cout << "Benchmarking sequential..." << std::flush;
        auto seq_result = benchmark_sequential(filename, size);
        seq_results.push_back(seq_result);
        std::cout << " " << static_cast<int>(seq_result.records_per_sec) << " rec/sec" << std::endl;
        
        std::cout << "Benchmarking parallel (4 threads)..." << std::flush;
        auto par_result = benchmark_parallel(filename, 4);
        par_results.push_back(par_result);
        std::cout << " " << static_cast<int>(par_result.records_per_sec) << " rec/sec" << std::endl;
        
        // Cleanup
        std::remove(filename.c_str());
        std::cout << std::endl;
    }
    
    // Print summary table
    std::cout << "=== Performance Summary ===" << std::endl;
    std::cout << "Records\t\tSeq (rec/s)\tPar (rec/s)\tSpeedup" << std::endl;
    std::cout << "-------\t\t-----------\t-----------\t-------" << std::endl;
    
    for (size_t i = 0; i < dataset_sizes.size(); ++i) {
        double speedup = par_results[i].records_per_sec / seq_results[i].records_per_sec;
        std::cout << dataset_sizes[i] << "\t\t"
                  << static_cast<int>(seq_results[i].records_per_sec) << "\t\t"
                  << static_cast<int>(par_results[i].records_per_sec) << "\t\t"
                  << speedup << "x" << std::endl;
    }
    
    // Check for O(n) scaling
    std::cout << "\n=== Scaling Analysis ===" << std::endl;
    if (seq_results.size() >= 2) {
        double first_rate = seq_results[0].records_per_sec;
        double last_rate = seq_results.back().records_per_sec;
        double variation = std::abs(last_rate - first_rate) / first_rate * 100;
        
        std::cout << "Sequential rate variation: " << variation << "%" << std::endl;
        if (variation < 20) {
            std::cout << "✓ Consistent O(n) scaling (< 20% variation)" << std::endl;
        } else {
            std::cout << "⚠ Performance degrades with scale (> 20% variation)" << std::endl;
        }
    }
}

// Test parallel scalability (thread count vs performance)
void test_parallel_scalability(const std::string& filename, const std::vector<size_t>& thread_counts) {
    std::cout << "\n=== Parallel Scalability Test ===" << std::endl;
    std::cout << "File: " << filename << std::endl;
    std::cout << "Thread counts: ";
    for (size_t count : thread_counts) {
        std::cout << count << " ";
    }
    std::cout << "\n" << std::endl;
    
    // Baseline: sequential
    std::cout << "Baseline (sequential)..." << std::flush;
    auto baseline = benchmark_sequential(filename, 0);
    std::cout << " " << static_cast<int>(baseline.records_per_sec) << " rec/sec" << std::endl;
    
    // Test with different thread counts
    std::cout << "\nThreads\tRate (rec/s)\tSpeedup\t\tEfficiency" << std::endl;
    std::cout << "-------\t------------\t-------\t\t----------" << std::endl;
    
    for (size_t threads : thread_counts) {
        auto result = benchmark_parallel(filename, threads);
        double speedup = result.records_per_sec / baseline.records_per_sec;
        double efficiency = speedup / threads * 100;
        
        std::cout << threads << "\t"
                  << static_cast<int>(result.records_per_sec) << "\t\t"
                  << speedup << "x\t\t"
                  << static_cast<int>(efficiency) << "%" << std::endl;
    }
}

int main(int argc, char* argv[]) {
    std::cout << "Boost.Genetics VCF Scalability Benchmark" << std::endl;
    std::cout << "=========================================" << std::endl;
    
    if (argc > 1 && std::string(argv[1]) == "--quick") {
        // Quick test with smaller datasets
        std::vector<size_t> sizes = {1000, 10000, 100000};
        test_scalability(sizes);
    } else if (argc > 1 && std::string(argv[1]) == "--parallel") {
        // Test parallel scaling with existing file
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " --parallel <vcf_file>" << std::endl;
            return 1;
        }
        std::string filename = argv[2];
        std::vector<size_t> thread_counts = {1, 2, 4, 6, 8};
        test_parallel_scalability(filename, thread_counts);
    } else {
        // Full scalability test
        std::cout << "\nRunning full scalability test (this may take several minutes)..." << std::endl;
        std::vector<size_t> sizes = {10000, 100000, 1000000};
        
        // Add 10M only if user explicitly requests
        if (argc > 1 && std::string(argv[1]) == "--full") {
            sizes.push_back(10000000);
        }
        
        test_scalability(sizes);
    }
    
    return 0;
}
