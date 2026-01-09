//  Memory profiling benchmark for Boost.Genetics VCF parser
//  Tracks peak memory usage and arena allocation efficiency

#include <boost/genetics/vcf.hpp>
#include <iostream>
#include <chrono>
#include <fstream>
#include <sys/resource.h>
#include <unistd.h>

using namespace boost::genetics;

// Get current memory usage in bytes (RSS)
size_t getCurrentRSS() {
#if defined(__APPLE__) && defined(__MACH__)
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    return static_cast<size_t>(rusage.ru_maxrss); // macOS returns bytes
#elif defined(__linux__)
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        if (line.substr(0, 6) == "VmRSS:") {
            size_t kb = std::stoul(line.substr(7));
            return kb * 1024; // Convert KB to bytes
        }
    }
    return 0;
#else
    return 0;
#endif
}

// Get peak memory usage in bytes
size_t getPeakRSS() {
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
    return static_cast<size_t>(rusage.ru_maxrss); // macOS returns bytes
#elif defined(__linux__)
    return static_cast<size_t>(rusage.ru_maxrss) * 1024; // Linux returns KB
#else
    return 0;
#endif
}

// Count number of records in a VCF file
size_t count_records(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) return 0;
    
    size_t count = 0;
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '#') {
            count++;
        }
    }
    return count;
}

// Benchmark memory usage for sequential reading
void benchmark_sequential_memory(const std::string& filename) {
    std::cout << "\n=== Sequential Reader Memory Benchmark ===" << std::endl;
    std::cout << "File: " << filename << std::endl;
    
    size_t record_count = count_records(filename);
    std::cout << "Records: " << record_count << std::endl;
    
    // Measure baseline memory
    size_t baseline_mem = getCurrentRSS();
    std::cout << "Baseline memory: " << (baseline_mem / 1024.0 / 1024.0) << " MB" << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    vcf::reader reader(filename);
    vcf::record rec;
    
    size_t records_read = 0;
    size_t peak_during_read = baseline_mem;
    
    while (reader.read_record(rec)) {
        records_read++;
        
        // Sample memory every 1000 records
        if (records_read % 1000 == 0) {
            size_t current = getCurrentRSS();
            if (current > peak_during_read) {
                peak_during_read = current;
            }
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    size_t final_mem = getCurrentRSS();
    size_t peak_mem = getPeakRSS();
    
    std::cout << "Records read: " << records_read << std::endl;
    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;
    std::cout << "Rate: " << (records_read / elapsed.count()) << " records/sec" << std::endl;
    std::cout << "\nMemory Usage:" << std::endl;
    std::cout << "  Peak RSS: " << (peak_mem / 1024.0 / 1024.0) << " MB" << std::endl;
    std::cout << "  Final RSS: " << (final_mem / 1024.0 / 1024.0) << " MB" << std::endl;
    std::cout << "  Memory overhead: " << ((final_mem - baseline_mem) / 1024.0 / 1024.0) << " MB" << std::endl;
    std::cout << "  Bytes per record: " << ((final_mem - baseline_mem) / static_cast<double>(records_read)) << std::endl;
}

// Benchmark memory usage for parallel reading
void benchmark_parallel_memory(const std::string& filename, size_t num_threads) {
    std::cout << "\n=== Parallel Reader Memory Benchmark ===" << std::endl;
    std::cout << "File: " << filename << std::endl;
    std::cout << "Threads: " << num_threads << std::endl;
    
    size_t record_count = count_records(filename);
    std::cout << "Records: " << record_count << std::endl;
    
    size_t baseline_mem = getCurrentRSS();
    std::cout << "Baseline memory: " << (baseline_mem / 1024.0 / 1024.0) << " MB" << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    vcf::parallel_reader reader(filename, num_threads);
    auto records = reader.read_all_unordered();
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    size_t final_mem = getCurrentRSS();
    size_t peak_mem = getPeakRSS();
    
    std::cout << "Records read: " << records.size() << std::endl;
    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;
    std::cout << "Rate: " << (records.size() / elapsed.count()) << " records/sec" << std::endl;
    std::cout << "\nMemory Usage:" << std::endl;
    std::cout << "  Peak RSS: " << (peak_mem / 1024.0 / 1024.0) << " MB" << std::endl;
    std::cout << "  Final RSS: " << (final_mem / 1024.0 / 1024.0) << " MB" << std::endl;
    std::cout << "  Memory overhead: " << ((final_mem - baseline_mem) / 1024.0 / 1024.0) << " MB" << std::endl;
    std::cout << "  Bytes per record: " << ((final_mem - baseline_mem) / static_cast<double>(records.size())) << std::endl;
    
    // Memory per thread estimation
    std::cout << "  Memory per thread: " << ((final_mem - baseline_mem) / num_threads / 1024.0 / 1024.0) << " MB" << std::endl;
}

// Benchmark arena allocator efficiency
void benchmark_arena_efficiency(const std::string& filename) {
    std::cout << "\n=== Arena Allocator Efficiency Benchmark ===" << std::endl;
    std::cout << "File: " << filename << std::endl;
    
    vcf::reader reader(filename);
    vcf::record rec;
    
    size_t records_read = 0;
    size_t total_string_bytes = 0;
    
    while (reader.read_record(rec)) {
        records_read++;
        
        // Estimate string data size per record
        total_string_bytes += rec.chrom().size();
        total_string_bytes += rec.id().size();
        total_string_bytes += rec.ref().size();
        for (const auto& alt : rec.alt()) {
            total_string_bytes += alt.size();
        }
        total_string_bytes += rec.filter().size();
    }
    
    std::cout << "Records analyzed: " << records_read << std::endl;
    std::cout << "Total string data: " << (total_string_bytes / 1024.0 / 1024.0) << " MB" << std::endl;
    std::cout << "Average string data per record: " << (total_string_bytes / static_cast<double>(records_read)) << " bytes" << std::endl;
    
    // Compare with actual memory usage
    size_t actual_mem = getCurrentRSS();
    std::cout << "\nArena efficiency:" << std::endl;
    std::cout << "  Theoretical minimum: " << (total_string_bytes / 1024.0 / 1024.0) << " MB" << std::endl;
    std::cout << "  Actual usage: " << (actual_mem / 1024.0 / 1024.0) << " MB" << std::endl;
    std::cout << "  Overhead factor: " << (static_cast<double>(actual_mem) / total_string_bytes) << "x" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <vcf_file> [num_threads]" << std::endl;
        return 1;
    }
    
    std::string filename = argv[1];
    size_t num_threads = (argc > 2) ? std::stoul(argv[2]) : 4;
    
    std::cout << "Boost.Genetics VCF Memory Profiling Benchmark" << std::endl;
    std::cout << "=============================================" << std::endl;
    
    // Run benchmarks
    benchmark_sequential_memory(filename);
    benchmark_arena_efficiency(filename);
    benchmark_parallel_memory(filename, num_threads);
    
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Peak memory: " << (getPeakRSS() / 1024.0 / 1024.0) << " MB" << std::endl;
    
    return 0;
}
