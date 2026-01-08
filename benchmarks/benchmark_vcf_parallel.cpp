#include <genetics/vcf.hpp>
#include <iostream>
#include <chrono>
#include <iomanip>

using namespace boost::genetics;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <vcf_file> [num_threads]\n";
        return 1;
    }
    
    std::string vcf_file = argv[1];
    std::size_t num_threads = (argc > 2) ? std::atoi(argv[2]) : 0; // 0 = auto-detect
    
    std::cout << "=== Parallel VCF Parsing Benchmark ===\n\n";
    std::cout << "Input file: " << vcf_file << "\n";
    std::cout << "Threads: " << (num_threads == 0 ? "auto (" + std::to_string(std::thread::hardware_concurrency()) + ")" : std::to_string(num_threads)) << "\n\n";
    
    // Get file size
    std::ifstream file_check(vcf_file, std::ios::binary | std::ios::ate);
    if (!file_check.is_open()) {
        std::cerr << "Error: Cannot open file " << vcf_file << "\n";
        return 1;
    }
    size_t file_size = file_check.tellg();
    file_check.close();
    
    // Sequential parsing (baseline)
    {
        std::cout << "Running sequential parser...\n";
        auto start = std::chrono::high_resolution_clock::now();
        
        vcf::reader reader(vcf_file);
        size_t count = 0;
        vcf::record rec;
        while (reader.read_record(rec)) {
            count++;
            // Access fields to ensure parsing
            (void)rec.chrom();
            (void)rec.pos();
            (void)rec.ref();
            (void)rec.alt().size();
            (void)rec.info();
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start).count();
        
        std::cout << "Sequential:                          : " << std::fixed << std::setprecision(3) 
                  << elapsed << "s, " << (int)(count / elapsed) << " rec/s, "
                  << std::setprecision(2) << (file_size / elapsed / (1024 * 1024)) << " MB/s\n";
        std::cout << "Records parsed: " << count << "\n\n";
    }
    
    // Parallel parsing
    {
        std::cout << "Running parallel parser...\n";
        auto start = std::chrono::high_resolution_clock::now();
        
        vcf::parallel_reader reader(vcf_file, num_threads);
        reader.read_header();
        std::vector<vcf::record> records = reader.read_all_unordered();
        
        // Access fields to ensure parsing
        for (auto& rec : records) {
            (void)rec.chrom();
            (void)rec.pos();
            (void)rec.ref();
            (void)rec.alt().size();
            (void)rec.info();
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start).count();
        
        std::cout << "Parallel (unordered):                : " << std::fixed << std::setprecision(3) 
                  << elapsed << "s, " << (int)(records.size() / elapsed) << " rec/s, "
                  << std::setprecision(2) << (file_size / elapsed / (1024 * 1024)) << " MB/s\n";
        std::cout << "Records parsed: " << records.size() << "\n";
    }
    
    return 0;
}
