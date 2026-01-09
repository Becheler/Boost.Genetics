//  VCF Writer performance benchmark
//  Measures write throughput and compares with reading

#include <boost/genetics/vcf.hpp>
#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>

using namespace boost::genetics;

// Benchmark sequential writing
void benchmark_sequential_write(const std::string& input_file, const std::string& output_file) {
    std::cout << "\n=== Sequential Write Benchmark ===" << std::endl;
    std::cout << "Input: " << input_file << std::endl;
    std::cout << "Output: " << output_file << std::endl;
    
    // First, read all records into memory
    std::cout << "Reading input file..." << std::endl;
    auto read_start = std::chrono::high_resolution_clock::now();
    
    vcf::reader reader(input_file);
    std::vector<vcf::record> records;
    vcf::record rec;
    
    while (reader.read_record(rec)) {
        records.push_back(rec);
    }
    
    auto read_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_time = read_end - read_start;
    
    std::cout << "Records read: " << records.size() << std::endl;
    std::cout << "Read time: " << read_time.count() << " seconds" << std::endl;
    std::cout << "Read rate: " << (records.size() / read_time.count()) << " records/sec" << std::endl;
    
    // Now benchmark writing
    std::cout << "\nWriting output file..." << std::endl;
    auto write_start = std::chrono::high_resolution_clock::now();
    
    vcf::writer writer(output_file);
    
    // Write header
    for (const auto& line : reader.header_lines()) {
        writer.add_header_line(line);
    }
    writer.write_header();
    
    // Write records
    for (const auto& record : records) {
        writer.write_record(record);
    }
    
    auto write_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> write_time = write_end - write_start;
    
    std::cout << "Records written: " << records.size() << std::endl;
    std::cout << "Write time: " << write_time.count() << " seconds" << std::endl;
    std::cout << "Write rate: " << (records.size() / write_time.count()) << " records/sec" << std::endl;
    
    // Compare read vs write performance
    std::cout << "\n=== Read vs Write Comparison ===" << std::endl;
    std::cout << "Read rate:  " << (records.size() / read_time.count()) << " records/sec" << std::endl;
    std::cout << "Write rate: " << (records.size() / write_time.count()) << " records/sec" << std::endl;
    std::cout << "Write/Read ratio: " << (write_time.count() / read_time.count()) << "x" << std::endl;
}

// Benchmark write throughput with generated data
void benchmark_synthetic_write(const std::string& output_file, size_t num_records) {
    std::cout << "\n=== Synthetic Data Write Benchmark ===" << std::endl;
    std::cout << "Output: " << output_file << std::endl;
    std::cout << "Records to generate: " << num_records << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    vcf::writer writer(output_file);
    
    // Write minimal header
    writer.add_header_line("##contig=<ID=chr1,length=248956422>");
    writer.add_header_line("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    writer.add_header_line("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    writer.write_header();
    
    // Generate and write records
    for (size_t i = 0; i < num_records; ++i) {
        std::vector<std::string> alt = {"G"};
        std::map<std::string, std::string> info;
        info["DP"] = std::to_string(50 + (i % 100));
        
        vcf::record rec("chr1", 10000 + i, ".", "A", alt, 30.0, "PASS", info);
        writer.write_record(rec);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;
    std::cout << "Rate: " << (num_records / elapsed.count()) << " records/sec" << std::endl;
    
    // Check output file size
    std::ifstream file(output_file, std::ios::ate | std::ios::binary);
    size_t file_size = file.tellg();
    std::cout << "Output file size: " << (file_size / 1024.0 / 1024.0) << " MB" << std::endl;
    std::cout << "Bytes per record: " << (file_size / static_cast<double>(num_records)) << std::endl;
    std::cout << "Write throughput: " << ((file_size / 1024.0 / 1024.0) / elapsed.count()) << " MB/sec" << std::endl;
}

// Benchmark round-trip (read + write)
void benchmark_roundtrip(const std::string& input_file, const std::string& output_file) {
    std::cout << "\n=== Round-trip Benchmark (Read + Write) ===" << std::endl;
    std::cout << "Input: " << input_file << std::endl;
    std::cout << "Output: " << output_file << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    vcf::reader reader(input_file);
    vcf::writer writer(output_file);
    
    // Copy header
    for (const auto& line : reader.header_lines()) {
        writer.add_header_line(line);
    }
    writer.write_header();
    
    // Stream records (no intermediate storage)
    vcf::record rec;
    size_t records_processed = 0;
    
    while (reader.read_record(rec)) {
        writer.write_record(rec);
        records_processed++;
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    std::cout << "Records processed: " << records_processed << std::endl;
    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;
    std::cout << "Throughput: " << (records_processed / elapsed.count()) << " records/sec" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_vcf_file> [output_vcf_file]" << std::endl;
        std::cerr << "   Or: " << argv[0] << " --synthetic <output_file> <num_records>" << std::endl;
        return 1;
    }
    
    std::cout << "Boost.Genetics VCF Writer Benchmark" << std::endl;
    std::cout << "====================================" << std::endl;
    
    if (std::string(argv[1]) == "--synthetic") {
        if (argc < 4) {
            std::cerr << "Error: --synthetic requires output file and record count" << std::endl;
            return 1;
        }
        std::string output_file = argv[2];
        size_t num_records = std::stoul(argv[3]);
        benchmark_synthetic_write(output_file, num_records);
    } else {
        std::string input_file = argv[1];
        std::string output_file = (argc > 2) ? argv[2] : "output_benchmark.vcf";
        
        benchmark_sequential_write(input_file, output_file);
        benchmark_roundtrip(input_file, "roundtrip_benchmark.vcf");
        
        // Cleanup
        std::remove(output_file.c_str());
        std::remove("roundtrip_benchmark.vcf");
    }
    
    return 0;
}
