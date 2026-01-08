//  Copyright (c) 2026 Boost.Genetics
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/genetics/vcf.hpp>
#include <iostream>
#include <vector>

int main() {
    using namespace boost::genetics;

    std::cout << "=== VCF Writer Example ===" << std::endl;
    
    // Create a VCF writer
    vcf::writer writer("example_output.vcf", "VCFv4.3");
    
    // Add header metadata
    writer.add_header_line("##reference=hg38");
    writer.add_contig("chr1", 248956422);
    writer.add_contig("chr2", 242193529);
    
    // Add INFO field definitions
    writer.add_info("DP", "1", "Integer", "Total Depth");
    writer.add_info("AF", "A", "Float", "Allele Frequency");
    writer.add_info("DB", "0", "Flag", "dbSNP membership");
    
    // Add FORMAT field definitions
    writer.add_format("GT", "1", "String", "Genotype");
    writer.add_format("GQ", "1", "Integer", "Genotype Quality");
    writer.add_format("DP", "1", "Integer", "Read Depth");
    
    // Set sample names
    std::vector<std::string> samples;
    samples.push_back("Sample1");
    samples.push_back("Sample2");
    writer.set_sample_names(samples);
    
    // Write header
    writer.write_header();
    
    // Create and write some variant records
    
    // SNP example
    {
        std::vector<std::string> alt;
        alt.push_back("G");
        
        std::map<std::string, std::string> info;
        info["DP"] = "100";
        info["AF"] = "0.5";
        info["DB"] = "";  // Flag field
        
        vcf::record rec("chr1", 12345, "rs123456", "A", alt, 30.0, "PASS", info);
        rec.set_format("GT:GQ:DP");
        
        // Add sample 1 data
        std::map<std::string, std::string> sample1;
        sample1["GT"] = "0/1";
        sample1["GQ"] = "99";
        sample1["DP"] = "45";
        rec.add_sample(sample1);
        
        // Add sample 2 data
        std::map<std::string, std::string> sample2;
        sample2["GT"] = "1/1";
        sample2["GQ"] = "80";
        sample2["DP"] = "55";
        rec.add_sample(sample2);
        
        writer.write_record(rec);
        std::cout << "Wrote SNP: chr1:12345 A>G" << std::endl;
    }
    
    // Insertion example
    {
        std::vector<std::string> alt;
        alt.push_back("TCGA");
        
        std::map<std::string, std::string> info;
        info["DP"] = "80";
        info["AF"] = "0.25";
        
        vcf::record rec("chr1", 67890, ".", "T", alt, 25.0, "PASS", info);
        rec.set_format("GT:GQ:DP");
        
        std::map<std::string, std::string> sample1;
        sample1["GT"] = "0/1";
        sample1["GQ"] = "70";
        sample1["DP"] = "40";
        rec.add_sample(sample1);
        
        std::map<std::string, std::string> sample2;
        sample2["GT"] = "0/0";
        sample2["GQ"] = "99";
        sample2["DP"] = "40";
        rec.add_sample(sample2);
        
        writer.write_record(rec);
        std::cout << "Wrote insertion: chr1:67890 T>TCGA" << std::endl;
    }
    
    // Deletion example
    {
        std::vector<std::string> alt;
        alt.push_back("C");
        
        std::map<std::string, std::string> info;
        info["DP"] = "120";
        info["AF"] = "0.75";
        
        vcf::record rec("chr2", 100000, "rs789012", "CGTA", alt, 40.0, "PASS", info);
        rec.set_format("GT:GQ:DP");
        
        std::map<std::string, std::string> sample1;
        sample1["GT"] = "1/1";
        sample1["GQ"] = "95";
        sample1["DP"] = "60";
        rec.add_sample(sample1);
        
        std::map<std::string, std::string> sample2;
        sample2["GT"] = "0/1";
        sample2["GQ"] = "85";
        sample2["DP"] = "60";
        rec.add_sample(sample2);
        
        writer.write_record(rec);
        std::cout << "Wrote deletion: chr2:100000 CGTA>C" << std::endl;
    }
    
    std::cout << "\n=== VCF Reader Example ===" << std::endl;
    
    // Read the file we just wrote
    try {
        vcf::reader reader("example_output.vcf");
        
        std::cout << "VCF Version: " << reader.version() << std::endl;
        std::cout << "Number of samples: " << reader.sample_names().size() << std::endl;
        std::cout << "Sample names: ";
        const std::vector<std::string>& sample_names = reader.sample_names();
        for (std::size_t i = 0; i < sample_names.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << sample_names[i];
        }
        std::cout << std::endl << std::endl;
        
        // Read all records
        std::vector<vcf::record> records = reader.read_all();
        std::cout << "Total variants: " << records.size() << std::endl << std::endl;
        
        // Display each variant
        for (std::size_t i = 0; i < records.size(); ++i) {
            const vcf::record& rec = records[i];
            std::cout << "Variant " << (i + 1) << ":" << std::endl;
            std::cout << "  Position: " << rec.chrom() << ":" << rec.pos() << std::endl;
            std::cout << "  ID: " << (rec.id().empty() ? "." : rec.id()) << std::endl;
            std::cout << "  REF: " << rec.ref() << std::endl;
            std::cout << "  ALT: ";
            const std::vector<std::string>& alts = rec.alt();
            for (std::size_t j = 0; j < alts.size(); ++j) {
                if (j > 0) std::cout << ",";
                std::cout << alts[j];
            }
            std::cout << std::endl;
            std::cout << "  QUAL: " << rec.qual() << std::endl;
            std::cout << "  FILTER: " << (rec.filter().empty() ? "." : rec.filter()) << std::endl;
            
            // Display INFO fields
            const std::map<std::string, std::string>& info = rec.info();
            std::cout << "  INFO: ";
            bool first = true;
            for (std::map<std::string, std::string>::const_iterator it = info.begin();
                 it != info.end(); ++it) {
                if (!first) std::cout << "; ";
                std::cout << it->first;
                if (!it->second.empty()) {
                    std::cout << "=" << it->second;
                }
                first = false;
            }
            std::cout << std::endl;
            
            // Display genotype data
            const auto& samples_data = rec.samples();
            const auto& format_fields = rec.get_format_fields();
            for (std::size_t j = 0; j < samples_data.size(); ++j) {
                std::cout << "  " << sample_names[j] << ": ";
                const auto& sample = samples_data[j];
                for (std::size_t k = 0; k < sample.size() && k < format_fields.size(); ++k) {
                    if (k > 0) std::cout << " ";
                    std::cout << format_fields[k] << "=" << sample[k];
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error reading VCF: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "VCF example completed successfully!" << std::endl;
    return 0;
}
