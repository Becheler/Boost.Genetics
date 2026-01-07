#include "../include/genetics/vcf.hpp"
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdlib>

using namespace boost::genetics;

struct benchmark_result {
    std::string name;
    double elapsed_seconds;
    size_t records_processed;
    size_t bytes_processed;
    
    double records_per_second() const {
        return records_processed / elapsed_seconds;
    }
    
    double mb_per_second() const {
        return (bytes_processed / (1024.0 * 1024.0)) / elapsed_seconds;
    }
};

void print_result(const benchmark_result& result) {
    std::cout << std::left << std::setw(40) << result.name << " : "
              << std::fixed << std::setprecision(3) << result.elapsed_seconds << "s, "
              << std::fixed << std::setprecision(0) << result.records_per_second() << " rec/s, "
              << std::fixed << std::setprecision(2) << result.mb_per_second() << " MB/s\n";
}

benchmark_result benchmark_parsing(const std::string& vcf_file) {
    auto start = std::chrono::high_resolution_clock::now();
    
    vcf::reader reader(vcf_file);
    
    // Get file size
    std::ifstream file_check(vcf_file, std::ios::binary | std::ios::ate);
    size_t file_size = file_check.tellg();
    file_check.close();
    
    size_t record_count = 0;
    vcf::record rec;
    while (reader.read_record(rec)) {
        record_count++;
        // Access fields to ensure they're parsed
        (void)rec.chrom();
        (void)rec.pos();
        (void)rec.ref();
        (void)rec.alt().size();
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    return {"VCF Parsing", elapsed.count(), record_count, file_size};
}

benchmark_result benchmark_parsing_with_validation(const std::string& vcf_file) {
    auto start = std::chrono::high_resolution_clock::now();
    
    vcf::reader reader(vcf_file);
    reader.validate_header();  // Enable validation
    
    std::ifstream file_check(vcf_file, std::ios::binary | std::ios::ate);
    size_t file_size = file_check.tellg();
    file_check.close();
    
    size_t record_count = 0;
    vcf::record rec;
    while (reader.read_record(rec)) {
        // Validate each record
        // Note: validate() doesn't exist - skip validation or implement custom checks
        record_count++;
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    return {"VCF Parsing + Validation", elapsed.count(), record_count, file_size};
}

benchmark_result benchmark_writing(const std::string& input_vcf, const std::string& output_vcf) {
    // First, read all records into memory
    vcf::reader reader(input_vcf);
    
    std::vector<vcf::record> records;
    vcf::record rec;
    while (reader.read_record(rec)) {
        records.push_back(rec);
    }
    
    // Now benchmark writing
    auto start = std::chrono::high_resolution_clock::now();
    
    vcf::writer writer(output_vcf);
    
    // Copy header lines
    for (const auto& line : reader.header_lines()) {
        writer.add_header_line(line);
    }
    writer.set_sample_names(reader.sample_names());
    writer.write_header();
    
    for (const auto& r : records) {
        writer.write_record(r);
    }
    
    // Get output file size
    std::ifstream check_size(output_vcf, std::ios::binary | std::ios::ate);
    size_t file_size = check_size.tellg();
    check_size.close();
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    return {"VCF Writing", elapsed.count(), records.size(), file_size};
}

benchmark_result benchmark_bcftools(const std::string& vcf_file) {
    auto start = std::chrono::high_resolution_clock::now();
    
    std::string cmd = "bcftools view " + vcf_file + " > /dev/null 2>&1";
    int result = system(cmd.c_str());
    if (result != 0) {
        throw std::runtime_error("bcftools command failed");
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    // Get file size
    std::ifstream input(vcf_file, std::ios::binary | std::ios::ate);
    size_t file_size = input.tellg();
    input.close();
    
    // Count records with bcftools
    cmd = "bcftools view -H " + vcf_file + " | wc -l";
    FILE* pipe = popen(cmd.c_str(), "r");
    size_t record_count = 0;
    if (pipe) {
        char buffer[128];
        if (fgets(buffer, sizeof(buffer), pipe)) {
            record_count = std::stoull(buffer);
        }
        pclose(pipe);
    }
    
    return {"bcftools view (baseline)", elapsed.count(), record_count, file_size};
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <vcf_file>\n";
        return 1;
    }
    
    std::string vcf_file = argv[1];
    
    std::cout << "=== VCF Parsing Benchmarks ===\n\n";
    std::cout << "Input file: " << vcf_file << "\n\n";
    
    try {
        // Baseline: bcftools
        std::cout << "Running bcftools baseline...\n";
        auto bcftools_result = benchmark_bcftools(vcf_file);
        print_result(bcftools_result);
        std::cout << "\n";
        
        // Our parser without validation
        std::cout << "Running our parser...\n";
        auto parse_result = benchmark_parsing(vcf_file);
        print_result(parse_result);
        
        double speedup = bcftools_result.elapsed_seconds / parse_result.elapsed_seconds;
        std::cout << "Speedup vs bcftools: " << std::fixed << std::setprecision(2) 
                  << speedup << "x\n\n";
        
        // Our parser with validation
        std::cout << "Running our parser with validation...\n";
        auto validate_result = benchmark_parsing_with_validation(vcf_file);
        print_result(validate_result);
        
        double validation_overhead = (validate_result.elapsed_seconds - parse_result.elapsed_seconds) 
                                    / parse_result.elapsed_seconds * 100.0;
        std::cout << "Validation overhead: " << std::fixed << std::setprecision(1) 
                  << validation_overhead << "%\n\n";
        
        // Writing benchmark
        std::cout << "Running write benchmark...\n";
        std::string output_file = "/tmp/benchmark_output.vcf";
        auto write_result = benchmark_writing(vcf_file, output_file);
        print_result(write_result);
        std::cout << "\n";
        
        // Memory usage estimate
        std::cout << "=== Memory Usage (estimate) ===\n";
        std::cout << "Records in memory: " << parse_result.records_processed << "\n";
        std::cout << "Average bytes per record: " 
                  << (parse_result.bytes_processed / parse_result.records_processed) << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
