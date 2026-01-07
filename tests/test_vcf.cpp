//  Copyright (c) 2026 Boost.Genetics
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"
#include <genetics/vcf.hpp>
#include <sstream>
#include <fstream>

using namespace boost::genetics;

TEST_CASE("VCF record creation and basic operations", "[vcf]") {
    SECTION("Default constructor") {
        vcf::record rec;
        REQUIRE(rec.chrom().empty());
        REQUIRE(rec.pos() == 0);
        REQUIRE(rec.id().empty());
        REQUIRE(rec.ref().empty());
        REQUIRE(rec.alt().empty());
        REQUIRE(rec.qual() == 0.0);
        REQUIRE(rec.filter().empty());
    }
    
    SECTION("Constructor with parameters") {
        std::vector<std::string> alt;
        alt.push_back("G");
        
        std::map<std::string, std::string> info;
        info["DP"] = "100";
        info["AF"] = "0.5";
        
        vcf::record rec("chr1", 12345, "rs123", "A", alt, 30.0, "PASS", info);
        
        REQUIRE(rec.chrom() == "chr1");
        REQUIRE(rec.pos() == 12345);
        REQUIRE(rec.id() == "rs123");
        REQUIRE(rec.ref() == "A");
        REQUIRE(rec.alt().size() == 1);
        REQUIRE(rec.alt()[0] == "G");
        REQUIRE(rec.qual() == 30.0);
        REQUIRE(rec.filter() == "PASS");
        REQUIRE(rec.info().size() == 2);
        REQUIRE(rec.info().find("DP")->second == "100");
        REQUIRE(rec.info().find("AF")->second == "0.5");
    }
}

TEST_CASE("VCF record to string conversion", "[vcf]") {
    SECTION("Basic SNP without samples") {
        std::vector<std::string> alt;
        alt.push_back("G");
        
        std::map<std::string, std::string> info;
        info["DP"] = "100";
        
        vcf::record rec("chr1", 12345, "rs123", "A", alt, 30.0, "PASS", info);
        
        std::string line = rec.to_string();
        REQUIRE(line.find("chr1") != std::string::npos);
        REQUIRE(line.find("12345") != std::string::npos);
        REQUIRE(line.find("rs123") != std::string::npos);
        REQUIRE(line.find("A") != std::string::npos);
        REQUIRE(line.find("G") != std::string::npos);
        REQUIRE(line.find("PASS") != std::string::npos);
        REQUIRE(line.find("DP=100") != std::string::npos);
    }
    
    SECTION("Multiple ALT alleles") {
        std::vector<std::string> alt;
        alt.push_back("G");
        alt.push_back("T");
        
        vcf::record rec("chr1", 100, ".", "A", alt, 0.0, ".", std::map<std::string, std::string>());
        
        std::string line = rec.to_string();
        REQUIRE(line.find("G,T") != std::string::npos);
    }
    
    SECTION("With sample data") {
        std::vector<std::string> alt;
        alt.push_back("G");
        
        vcf::record rec("chr1", 12345, ".", "A", alt, 30.0, "PASS", std::map<std::string, std::string>());
        rec.set_format("GT:DP:GQ");
        
        std::map<std::string, std::string> sample1;
        sample1["GT"] = "0/1";
        sample1["DP"] = "50";
        sample1["GQ"] = "99";
        rec.add_sample(sample1);
        
        std::vector<std::string> sample_names;
        sample_names.push_back("Sample1");
        std::string line = rec.to_string(sample_names);
        
        REQUIRE(line.find("GT:DP:GQ") != std::string::npos);
        REQUIRE(line.find("0/1:50:99") != std::string::npos);
    }
}

TEST_CASE("VCF writer functionality", "[vcf]") {
    const std::string test_file = "test_output.vcf";
    
    SECTION("Write header and basic record") {
        {
            vcf::writer writer(test_file, "VCFv4.3");
            writer.add_contig("chr1", 248956422);
            writer.add_info("DP", "1", "Integer", "Total Depth");
            writer.add_format("GT", "1", "String", "Genotype");
            
            std::vector<std::string> samples;
            samples.push_back("Sample1");
            writer.set_sample_names(samples);
            
            writer.write_header();
            
            std::vector<std::string> alt;
            alt.push_back("G");
            
            std::map<std::string, std::string> info;
            info["DP"] = "100";
            
            vcf::record rec("chr1", 12345, "rs123", "A", alt, 30.0, "PASS", info);
            rec.set_format("GT");
            
            std::map<std::string, std::string> sample_data;
            sample_data["GT"] = "0/1";
            rec.add_sample(sample_data);
            
            writer.write_record(rec);
        }
        
        // Verify file was created and has content
        std::ifstream check(test_file.c_str());
        REQUIRE(check.is_open());
        
        std::string first_line;
        std::getline(check, first_line);
        REQUIRE(first_line == "##fileformat=VCFv4.3");
        
        check.close();
        std::remove(test_file.c_str());
    }
}

TEST_CASE("VCF reader functionality", "[vcf]") {
    const std::string test_file = "test_read.vcf";
    
    // Create a test VCF file
    {
        std::ofstream out(test_file.c_str());
        out << "##fileformat=VCFv4.3\n";
        out << "##contig=<ID=chr1,length=248956422>\n";
        out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
        out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n";
        out << "chr1\t12345\trs123\tA\tG\t30\tPASS\tDP=100\tGT\t0/1\n";
        out << "chr1\t67890\t.\tT\tC\t40\tPASS\tDP=80\tGT\t1/1\n";
        out.close();
    }
    
    SECTION("Read header information") {
        vcf::reader reader(test_file);
        
        REQUIRE(reader.version() == "VCFv4.3");
        REQUIRE(reader.sample_names().size() == 1);
        REQUIRE(reader.sample_names()[0] == "Sample1");
        REQUIRE(reader.header_lines().size() > 0);
    }
    
    SECTION("Read records") {
        vcf::reader reader(test_file);
        
        vcf::record rec1;
        REQUIRE(reader.read_record(rec1));
        REQUIRE(rec1.chrom() == "chr1");
        REQUIRE(rec1.pos() == 12345);
        REQUIRE(rec1.id() == "rs123");
        REQUIRE(rec1.ref() == "A");
        REQUIRE(rec1.alt().size() == 1);
        REQUIRE(rec1.alt()[0] == "G");
        REQUIRE(rec1.qual() == 30.0);
        REQUIRE(rec1.filter() == "PASS");
        
        vcf::record rec2;
        REQUIRE(reader.read_record(rec2));
        REQUIRE(rec2.chrom() == "chr1");
        REQUIRE(rec2.pos() == 67890);
        REQUIRE(rec2.id().empty());
        
        vcf::record rec3;
        REQUIRE_FALSE(reader.read_record(rec3)); // No more records
    }
    
    SECTION("Read all records") {
        vcf::reader reader(test_file);
        std::vector<vcf::record> records = reader.read_all();
        
        REQUIRE(records.size() == 2);
        REQUIRE(records[0].pos() == 12345);
        REQUIRE(records[1].pos() == 67890);
    }
    
    std::remove(test_file.c_str());
}

TEST_CASE("VCF record INFO field operations", "[vcf]") {
    vcf::record rec;
    
    SECTION("Add INFO fields") {
        rec.add_info("DP", "100");
        rec.add_info("AF", "0.5");
        rec.add_info("DB", ""); // Flag field
        
        const std::map<std::string, std::string>& info = rec.info();
        REQUIRE(info.size() == 3);
        REQUIRE(info.find("DP")->second == "100");
        REQUIRE(info.find("AF")->second == "0.5");
        REQUIRE(info.find("DB")->second.empty());
    }
}

TEST_CASE("VCF variant type detection", "[vcf]") {
    SECTION("SNP") {
        std::vector<std::string> alt;
        alt.push_back("G");
        vcf::record rec("chr1", 100, ".", "A", alt, 0.0, ".", std::map<std::string, std::string>());
        
        REQUIRE(rec.ref().length() == 1);
        REQUIRE(rec.alt()[0].length() == 1);
    }
    
    SECTION("Insertion") {
        std::vector<std::string> alt;
        alt.push_back("AGCT");
        vcf::record rec("chr1", 100, ".", "A", alt, 0.0, ".", std::map<std::string, std::string>());
        
        REQUIRE(rec.ref().length() < rec.alt()[0].length());
    }
    
    SECTION("Deletion") {
        std::vector<std::string> alt;
        alt.push_back("A");
        vcf::record rec("chr1", 100, ".", "AGCT", alt, 0.0, ".", std::map<std::string, std::string>());
        
        REQUIRE(rec.ref().length() > rec.alt()[0].length());
    }
}

TEST_CASE("VCF sample data operations", "[vcf]") {
    vcf::record rec;
    rec.set_format("GT:DP:GQ");
    
    SECTION("Add multiple samples") {
        std::map<std::string, std::string> sample1;
        sample1["GT"] = "0/1";
        sample1["DP"] = "50";
        sample1["GQ"] = "99";
        rec.add_sample(sample1);
        
        std::map<std::string, std::string> sample2;
        sample2["GT"] = "1/1";
        sample2["DP"] = "60";
        sample2["GQ"] = "80";
        rec.add_sample(sample2);
        
        REQUIRE(rec.samples().size() == 2);
        REQUIRE(rec.samples()[0].find("GT")->second == "0/1");
        REQUIRE(rec.samples()[1].find("GT")->second == "1/1");
    }
}

TEST_CASE("VCF error handling with line numbers", "[vcf]") {
    const std::string test_file = "test_error.vcf";
    
    SECTION("Invalid VCF record throws vcf_parse_error with line number") {
        // Create a malformed VCF file
        {
            std::ofstream out(test_file.c_str());
            out << "##fileformat=VCFv4.3\n";
            out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
            out << "chr1\t12345\trs123\tA\n";  // Missing columns
            out.close();
        }
        
        try {
            vcf::reader reader(test_file);
            vcf::record rec;
            reader.read_record(rec);
            REQUIRE(false); // Should not reach here
        } catch (const vcf::vcf_parse_error& e) {
            REQUIRE(e.line() == 3); // Error on line 3
            REQUIRE(std::string(e.what()).find("Line 3") != std::string::npos);
        }
        
        std::remove(test_file.c_str());
    }
}

TEST_CASE("VCF missing genotype values (Spec Example 2.5)", "[vcf]") {
    const std::string test_file = "test_missing_gt.vcf";
    
    // Create test file matching spec example
    {
        std::ofstream out(test_file.c_str());
        out << "##fileformat=VCFv4.3\n";
        out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        out << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
        out << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";
        out << "chr1\t500\t.\tA\tG\t30\tPASS\t.\tGT:GQ:DP\t0/1:.:10\t./.:20:.\n";
        out.close();
    }
    
    vcf::reader reader(test_file);
    vcf::record rec;
    REQUIRE(reader.read_record(rec));
    
    REQUIRE(rec.chrom() == "chr1");
    REQUIRE(rec.pos() == 500);
    REQUIRE(rec.samples().size() == 2);
    
    // S1: GT=0/1, GQ=., DP=10
    REQUIRE(rec.samples()[0].find("GT")->second == "0/1");
    REQUIRE(rec.samples()[0].find("GQ")->second == ".");
    REQUIRE(rec.samples()[0].find("DP")->second == "10");
    
    // S2: GT=./., GQ=20, DP=.
    REQUIRE(rec.samples()[1].find("GT")->second == "./.");
    REQUIRE(rec.samples()[1].find("GQ")->second == "20");
    REQUIRE(rec.samples()[1].find("DP")->second == ".");
    
    std::remove(test_file.c_str());
}

TEST_CASE("VCF phased genotypes (Spec Example 2.6)", "[vcf]") {
    const std::string test_file = "test_phased.vcf";
    
    {
        std::ofstream out(test_file.c_str());
        out << "##fileformat=VCFv4.3\n";
        out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        out << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n";
        out << "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT:PS\t0|1:100\n";
        out << "chr1\t200\t.\tC\tT\t30\tPASS\t.\tGT:PS\t1|0:100\n";
        out << "chr1\t300\t.\tG\tA\t30\tPASS\t.\tGT\t0/1\n";
        out.close();
    }
    
    vcf::reader reader(test_file);
    
    vcf::record rec1;
    REQUIRE(reader.read_record(rec1));
    REQUIRE(rec1.samples()[0].find("GT")->second == "0|1");
    REQUIRE(rec1.samples()[0].find("PS")->second == "100");
    
    vcf::record rec2;
    REQUIRE(reader.read_record(rec2));
    REQUIRE(rec2.samples()[0].find("GT")->second == "1|0");
    REQUIRE(rec2.samples()[0].find("PS")->second == "100");
    
    vcf::record rec3;
    REQUIRE(reader.read_record(rec3));
    REQUIRE(rec3.samples()[0].find("GT")->second == "0/1");
    // No PS field in this record
    
    std::remove(test_file.c_str());
}

TEST_CASE("VCF multiple ALT alleles", "[vcf]") {
    const std::string test_file = "test_multi_alt.vcf";
    
    {
        std::ofstream out(test_file.c_str());
        out << "##fileformat=VCFv4.3\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        out << "chr1\t600\t.\tA\tG,T\t50\tPASS\t.\n";
        out.close();
    }
    
    vcf::reader reader(test_file);
    vcf::record rec;
    REQUIRE(reader.read_record(rec));
    
    REQUIRE(rec.alt().size() == 2);
    REQUIRE(rec.alt()[0] == "G");
    REQUIRE(rec.alt()[1] == "T");
    
    std::remove(test_file.c_str());
}

TEST_CASE("VCF Phase 2: Percent encoding/decoding", "[vcf][phase2]") {
    SECTION("Encode special characters") {
        std::string original = "value with spaces";
        std::string encoded = vcf::percent_encoding::encode(original);
        REQUIRE(encoded.find("%20") != std::string::npos);
        
        std::string decoded = vcf::percent_encoding::decode(encoded);
        REQUIRE(decoded == original);
    }
    
    SECTION("Encode semicolon and equals") {
        std::string original = "key=value;other";
        std::string encoded = vcf::percent_encoding::encode(original);
        REQUIRE(encoded.find("%3D") != std::string::npos); // =
        REQUIRE(encoded.find("%3B") != std::string::npos); // ;
        
        std::string decoded = vcf::percent_encoding::decode(encoded);
        REQUIRE(decoded == original);
    }
}

TEST_CASE("VCF Phase 2: Parse structured headers", "[vcf][phase2]") {
    SECTION("Parse INFO header") {
        std::string content = "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">";
        std::map<std::string, std::string> fields = vcf::parse_structured_header(content);
        
        REQUIRE(fields["ID"] == "DP");
        REQUIRE(fields["Number"] == "1");
        REQUIRE(fields["Type"] == "Integer");
        REQUIRE(fields["Description"] == "Total Depth");
    }
    
    SECTION("Parse contig header with length") {
        std::string content = "<ID=chr1,length=248956422>";
        std::map<std::string, std::string> fields = vcf::parse_structured_header(content);
        
        REQUIRE(fields["ID"] == "chr1");
        REQUIRE(fields["length"] == "248956422");
    }
}

TEST_CASE("VCF Phase 2: Reader parses metadata", "[vcf][phase2]") {
    const std::string test_file = "test_metadata.vcf";
    
    {
        std::ofstream out(test_file.c_str());
        out << "##fileformat=VCFv4.3\n";
        out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
        out << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n";
        out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        out << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
        out << "##FILTER=<ID=q10,Description=\"Quality below 10\">\n";
        out << "##contig=<ID=chr1,length=248956422>\n";
        out << "##contig=<ID=chr2,length=242193529>\n";
        out << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        out << "chr1\t100\t.\tA\tG\t30\tPASS\tDP=50;AF=0.5\n";
        out.close();
    }
    
    vcf::reader reader(test_file);
    
    SECTION("INFO metadata parsed") {
        const std::map<std::string, vcf::info_meta>& info_meta = reader.info_metadata();
        REQUIRE(info_meta.size() == 2);
        
        REQUIRE(info_meta.count("DP") == 1);
        REQUIRE(info_meta.at("DP").id == "DP");
        REQUIRE(info_meta.at("DP").number == "1");
        REQUIRE(info_meta.at("DP").type == "Integer");
        REQUIRE(info_meta.at("DP").description == "Total Depth");
        
        REQUIRE(info_meta.count("AF") == 1);
        REQUIRE(info_meta.at("AF").number == "A");
        REQUIRE(info_meta.at("AF").type == "Float");
    }
    
    SECTION("FORMAT metadata parsed") {
        const std::map<std::string, vcf::format_meta>& format_meta = reader.format_metadata();
        REQUIRE(format_meta.size() == 2);
        
        REQUIRE(format_meta.count("GT") == 1);
        REQUIRE(format_meta.at("GT").id == "GT");
        REQUIRE(format_meta.at("GT").type == "String");
        
        REQUIRE(format_meta.count("GQ") == 1);
    }
    
    SECTION("FILTER metadata parsed") {
        const std::map<std::string, vcf::filter_meta>& filter_meta = reader.filter_metadata();
        REQUIRE(filter_meta.size() == 1);
        
        REQUIRE(filter_meta.count("q10") == 1);
        REQUIRE(filter_meta.at("q10").description == "Quality below 10");
    }
    
    SECTION("Contig metadata parsed") {
        const std::map<std::string, vcf::contig_meta>& contig_meta = reader.contig_metadata();
        REQUIRE(contig_meta.size() == 2);
        
        REQUIRE(contig_meta.count("chr1") == 1);
        REQUIRE(contig_meta.at("chr1").length == 248956422);
        
        REQUIRE(contig_meta.count("chr2") == 1);
        REQUIRE(contig_meta.at("chr2").length == 242193529);
    }
    
    SECTION("ALT metadata parsed") {
        const std::map<std::string, vcf::alt_meta>& alt_meta = reader.alt_metadata();
        REQUIRE(alt_meta.size() == 1);
        
        REQUIRE(alt_meta.count("DEL") == 1);
        REQUIRE(alt_meta.at("DEL").description == "Deletion");
    }
    
    std::remove(test_file.c_str());
}

TEST_CASE("VCF Phase 2: Record validation", "[vcf][phase2]") {
    const std::string test_file = "test_validate.vcf";
    
    {
        std::ofstream out(test_file.c_str());
        out << "##fileformat=VCFv4.3\n";
        out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        out << "chr1\t100\t.\tA\tG\t30\tPASS\tDP=50\n";
        out.close();
    }
    
    vcf::reader reader(test_file);
    vcf::record rec;
    REQUIRE(reader.read_record(rec));
    
    SECTION("Valid INFO validates successfully") {
        REQUIRE(rec.validate_info(reader.info_metadata()));
    }
    
    SECTION("Invalid INFO fails validation") {
        rec.add_info("UNDEFINED_KEY", "value");
        REQUIRE_FALSE(rec.validate_info(reader.info_metadata()));
    }
    
    SECTION("Variant type detection") {
        REQUIRE(rec.is_snp());
        REQUIRE_FALSE(rec.is_insertion());
        REQUIRE_FALSE(rec.is_deletion());
    }
    
    std::remove(test_file.c_str());
}

// ============================================================================
// PHASE 3: STRUCTURAL VARIANT TESTS
// ============================================================================

TEST_CASE("Symbolic allele detection", "[vcf][phase3]") {
    SECTION("DEL symbolic allele") {
        REQUIRE(vcf::is_symbolic_allele("<DEL>"));
        REQUIRE(vcf::get_symbolic_id("<DEL>") == "DEL");
        REQUIRE(vcf::string_to_sv_type("DEL") == vcf::sv_type::DEL);
    }
    
    SECTION("INS symbolic allele") {
        REQUIRE(vcf::is_symbolic_allele("<INS>"));
        REQUIRE(vcf::get_symbolic_id("<INS>") == "INS");
    }
    
    SECTION("DUP symbolic allele") {
        REQUIRE(vcf::is_symbolic_allele("<DUP>"));
        REQUIRE(vcf::get_symbolic_id("<DUP>") == "DUP");
    }
    
    SECTION("INV symbolic allele") {
        REQUIRE(vcf::is_symbolic_allele("<INV>"));
        REQUIRE(vcf::get_symbolic_id("<INV>") == "INV");
    }
    
    SECTION("CNV symbolic allele") {
        REQUIRE(vcf::is_symbolic_allele("<CNV>"));
        REQUIRE(vcf::get_symbolic_id("<CNV>") == "CNV");
    }
    
    SECTION("Non-symbolic allele") {
        REQUIRE_FALSE(vcf::is_symbolic_allele("A"));
        REQUIRE_FALSE(vcf::is_symbolic_allele("ATG"));
        REQUIRE(vcf::get_symbolic_id("A").empty());
    }
}

TEST_CASE("Simple deletion with symbolic allele", "[vcf][phase3]") {
    const std::string test_file = "test_sv_deletion.vcf";
    {
        std::ofstream out(test_file.c_str());
        out << "##fileformat=VCFv4.3\n";
        out << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
        out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV\">\n";
        out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n";
        out << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length difference\">\n";
        out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n";
        out << "chr1\t1000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2000;SVLEN=-1000\tGT\t1/1\n";
        out.close();
    }
    
    vcf::reader reader(test_file);
    vcf::record rec;
    REQUIRE(reader.read_record(rec));
    
    SECTION("Detect symbolic SV") {
        REQUIRE(rec.is_symbolic_sv());
        REQUIRE_FALSE(rec.is_breakend());
    }
    
    SECTION("SV type detection") {
        REQUIRE(rec.get_sv_type() == vcf::sv_type::DEL);
    }
    
    SECTION("SV info parsing") {
        vcf::sv_info sv = rec.get_sv_info();
        REQUIRE(sv.type == vcf::sv_type::DEL);
        REQUIRE(sv.end == 2000);
        REQUIRE(sv.svlen == -1000);
        REQUIRE_FALSE(sv.imprecise);
    }
    
    SECTION("Helper methods") {
        REQUIRE(rec.get_end_pos() == 2000);
        REQUIRE(rec.get_svlen() == -1000);
        REQUIRE_FALSE(rec.is_imprecise());
    }
    
    std::remove(test_file.c_str());
}

TEST_CASE("Imprecise deletion", "[vcf][phase3]") {
    const std::string test_file = "test_sv_imprecise.vcf";
    {
        std::ofstream out(test_file.c_str());
        out << "##fileformat=VCFv4.3\n";
        out << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
        out << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise variant\">\n";
        out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV\">\n";
        out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n";
        out << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length difference\">\n";
        out << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"CI around POS\">\n";
        out << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"CI around END\">\n";
        out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n";
        out << "chr1\t5000\t.\tN\t<DEL>\t6\tPASS\tIMPRECISE;SVTYPE=DEL;END=5205;SVLEN=-205;CIPOS=-56,20;CIEND=-10,62\tGT\t0/1\n";
        out.close();
    }
    
    vcf::reader reader(test_file);
    vcf::record rec;
    REQUIRE(reader.read_record(rec));
    
    SECTION("Imprecise flag detection") {
        vcf::sv_info sv = rec.get_sv_info();
        REQUIRE(sv.imprecise);
        REQUIRE(rec.is_imprecise());
    }
    
    SECTION("Confidence intervals") {
        vcf::sv_info sv = rec.get_sv_info();
        REQUIRE(sv.cipos.first == -56);
        REQUIRE(sv.cipos.second == 20);
        REQUIRE(sv.ciend.first == -10);
        REQUIRE(sv.ciend.second == 62);
    }
    
    std::remove(test_file.c_str());
}

TEST_CASE("Breakend notation parsing", "[vcf][phase3]") {
    SECTION("Pattern: t[p[ - forward strand, joined after") {
        vcf::breakend bnd;
        REQUIRE(vcf::parse_breakend("G[chr17:198982[", bnd));
        REQUIRE(bnd.novel_sequence == "G");
        REQUIRE(bnd.mate_chr == "chr17");
        REQUIRE(bnd.mate_pos == 198982);
        REQUIRE(bnd.orientation == vcf::breakend_orientation::FORWARD);
        REQUIRE(bnd.position == vcf::breakend_position::AFTER);
    }
    
    SECTION("Pattern: t]p] - reverse complement, joined after") {
        vcf::breakend bnd;
        REQUIRE(vcf::parse_breakend("T]chr13:123456]", bnd));
        REQUIRE(bnd.novel_sequence == "T");
        REQUIRE(bnd.mate_chr == "chr13");
        REQUIRE(bnd.mate_pos == 123456);
        REQUIRE(bnd.orientation == vcf::breakend_orientation::REVERSE);
        REQUIRE(bnd.position == vcf::breakend_position::AFTER);
    }
    
    SECTION("Pattern: ]p]t - forward strand, joined before") {
        vcf::breakend bnd;
        REQUIRE(vcf::parse_breakend("]chr2:321682]T", bnd));
        REQUIRE(bnd.novel_sequence == "T");
        REQUIRE(bnd.mate_chr == "chr2");
        REQUIRE(bnd.mate_pos == 321682);
        REQUIRE(bnd.orientation == vcf::breakend_orientation::FORWARD);
        REQUIRE(bnd.position == vcf::breakend_position::BEFORE);
    }
    
    SECTION("Pattern: [p[t - reverse complement, joined before") {
        vcf::breakend bnd;
        REQUIRE(vcf::parse_breakend("[chr13:123456[C", bnd));
        REQUIRE(bnd.novel_sequence == "C");
        REQUIRE(bnd.mate_chr == "chr13");
        REQUIRE(bnd.mate_pos == 123456);
        REQUIRE(bnd.orientation == vcf::breakend_orientation::REVERSE);
        REQUIRE(bnd.position == vcf::breakend_position::BEFORE);
    }
    
    SECTION("Non-breakend notation") {
        vcf::breakend bnd;
        REQUIRE_FALSE(vcf::parse_breakend("A", bnd));
        REQUIRE_FALSE(vcf::parse_breakend("<DEL>", bnd));
        REQUIRE_FALSE(vcf::parse_breakend("ATG", bnd));
    }
}

TEST_CASE("Breakend VCF record", "[vcf][phase3]") {
    const std::string test_file = "test_breakend.vcf";
    {
        std::ofstream out(test_file.c_str());
        out << "##fileformat=VCFv4.3\n";
        out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV\">\n";
        out << "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"Mate breakend ID\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        out << "chr2\t321682\tbnd_V\tT\t]chr13:123456]T\t6\tPASS\tSVTYPE=BND;MATEID=bnd_U\n";
        out << "chr13\t123456\tbnd_U\tC\tC[chr2:321682[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_V\n";
        out.close();
    }
    
    vcf::reader reader(test_file);
    vcf::record rec;
    REQUIRE(reader.read_record(rec));
    
    SECTION("Detect breakend notation") {
        REQUIRE(rec.is_breakend());
        REQUIRE_FALSE(rec.is_symbolic_sv());
    }
    
    SECTION("SV type is BND") {
        REQUIRE(rec.get_sv_type() == vcf::sv_type::BND);
    }
    
    SECTION("Parse breakend") {
        vcf::breakend bnd;
        REQUIRE(rec.get_breakend(bnd));
        REQUIRE(bnd.novel_sequence == "T");
        REQUIRE(bnd.mate_chr == "chr13");
        REQUIRE(bnd.mate_pos == 123456);
        REQUIRE(bnd.orientation == vcf::breakend_orientation::FORWARD);
        REQUIRE(bnd.position == vcf::breakend_position::BEFORE);
    }
    
    SECTION("MATEID parsing") {
        vcf::sv_info sv = rec.get_sv_info();
        REQUIRE(sv.mateid.size() == 1);
        REQUIRE(sv.mateid[0] == "bnd_U");
    }
    
    std::remove(test_file.c_str());
}

TEST_CASE("Multiple SV types", "[vcf][phase3]") {
    SECTION("Duplication") {
        REQUIRE(vcf::string_to_sv_type("DUP") == vcf::sv_type::DUP);
        REQUIRE(vcf::sv_type_to_string(vcf::sv_type::DUP) == "DUP");
    }
    
    SECTION("Inversion") {
        REQUIRE(vcf::string_to_sv_type("INV") == vcf::sv_type::INV);
        REQUIRE(vcf::sv_type_to_string(vcf::sv_type::INV) == "INV");
    }
    
    SECTION("Copy number variation") {
        REQUIRE(vcf::string_to_sv_type("CNV") == vcf::sv_type::CNV);
        REQUIRE(vcf::sv_type_to_string(vcf::sv_type::CNV) == "CNV");
    }
    
    SECTION("Unknown type") {
        REQUIRE(vcf::string_to_sv_type("UNKNOWN") == vcf::sv_type::UNKNOWN);
        REQUIRE(vcf::sv_type_to_string(vcf::sv_type::UNKNOWN) == "UNKNOWN");
    }
}
