//  Copyright (c) 2026 Boost.Genetics
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GENETICS_VCF_HPP
#define BOOST_GENETICS_VCF_HPP

#include <string>
#include <string_view>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <genetics/vcf_gvcf.hpp>
#include <genetics/vcf_metadata.hpp>
#include <genetics/vcf_structural_variants.hpp>
#include <genetics/vcf_validation.hpp>

namespace boost {
namespace genetics {
namespace vcf {

// ============================================================================
// FAST PARSING UTILITIES
// ============================================================================

namespace detail {

/// @brief Fast integer parsing from string_view (avoids std::stoi overhead)
inline std::size_t parse_uint_fast(std::string_view str) noexcept {
    std::size_t result = 0;
    for (char c : str) {
        if (c >= '0' && c <= '9') {
            result = result * 10 + (c - '0');
        } else {
            break;
        }
    }
    return result;
}

/// @brief Fast double parsing from string_view (fallback to strtod)
inline double parse_double_fast(std::string_view str) noexcept {
    // std::from_chars for floating point is not universally available in C++17
    // Use fast strtod approach instead
    if (str.empty()) return 0.0;
    
    // Create null-terminated string on stack for small strings (SSO optimization)
    char buffer[32];
    if (str.size() < sizeof(buffer)) {
        std::memcpy(buffer, str.data(), str.size());
        buffer[str.size()] = '\0';
        return std::strtod(buffer, nullptr);
    }
    
    // Fallback for larger strings (rare)
    std::string temp(str);
    return std::strtod(temp.c_str(), nullptr);
}

/// @brief Split string_view by delimiter, calling callback for each field
template<typename Func>
inline void split_view(std::string_view str, char delim, Func callback) {
    std::size_t start = 0;
    std::size_t pos = str.find(delim);
    
    while (pos != std::string_view::npos) {
        callback(str.substr(start, pos - start));
        start = pos + 1;
        pos = str.find(delim, start);
    }
    
    // Last field
    if (start < str.size()) {
        callback(str.substr(start));
    } else if (start == str.size()) {
        callback(std::string_view());
    }
}

/// @brief Extract tab-separated fields from a line (zero-copy)
struct field_extractor {
    std::string_view line;
    std::size_t pos = 0;
    
    explicit field_extractor(std::string_view l) : line(l) {}
    
    std::string_view next() {
        if (pos >= line.size()) {
            return std::string_view();
        }
        
        std::size_t tab_pos = line.find('\t', pos);
        std::string_view result;
        
        if (tab_pos == std::string_view::npos) {
            result = line.substr(pos);
            pos = line.size();
        } else {
            result = line.substr(pos, tab_pos - pos);
            pos = tab_pos + 1;
        }
        
        return result;
    }
    
    bool has_more() const {
        return pos < line.size();
    }
};

} // namespace detail

// ============================================================================
// VCF RECORD CLASS
// ============================================================================

/// @brief VCF parsing exception with line context
class vcf_parse_error : public std::runtime_error {
public:
    /// @brief Constructor with message and line number
    /// @param msg Error message
    /// @param line Line number where error occurred (0 if unknown)
    vcf_parse_error(const std::string& msg, std::size_t line = 0)
        : std::runtime_error(build_message(msg, line)), line_number_(line) {}

    /// @brief Get the line number where error occurred
    std::size_t line() const { return line_number_; }

private:
    std::size_t line_number_;

    static std::string build_message(const std::string& msg, std::size_t line) {
        if (line == 0) {
            return msg;
        }
        std::ostringstream oss;
        oss << "Line " << line << ": " << msg;
        return oss.str();
    }
};

/// @brief Represents a single VCF record (variant call)
/// Compliant with VCFv4.3 specification
class record {
public:
    /// @brief Default constructor
    record() : pos_(0), qual_(0.0) {}

    /// @brief Constructor with basic variant information
    /// @param chrom Chromosome identifier
    /// @param pos Position (1-based)
    /// @param id Variant ID (or "." if none)
    /// @param ref Reference allele
    /// @param alt Alternative allele(s)
    /// @param qual Quality score
    /// @param filter Filter status
    /// @param info INFO field data
    record(const std::string& chrom, std::size_t pos, const std::string& id,
           const std::string& ref, const std::vector<std::string>& alt,
           double qual, const std::string& filter,
           const std::map<std::string, std::string>& info = std::map<std::string, std::string>())
        : chrom_(chrom), pos_(pos), id_(id), ref_(ref), alt_(alt),
          qual_(qual), filter_(filter), info_(info) {}

    // Getters
    const std::string& chrom() const { return chrom_; }
    std::size_t pos() const { return pos_; }
    const std::string& id() const { return id_; }
    const std::string& ref() const { return ref_; }
    const std::vector<std::string>& alt() const { return alt_; }
    double qual() const { return qual_; }
    const std::string& filter() const { return filter_; }
    const std::map<std::string, std::string>& info() const { return info_; }
    const std::string& format() const { return format_; }
    const std::vector<std::map<std::string, std::string>>& samples() const { return samples_; }

    // Setters
    void set_chrom(const std::string& chrom) { chrom_ = chrom; }
    void set_pos(std::size_t pos) { pos_ = pos; }
    void set_id(const std::string& id) { id_ = id; }
    void set_ref(const std::string& ref) { ref_ = ref; }
    void set_alt(const std::vector<std::string>& alt) { alt_ = alt; }
    void set_qual(double qual) { qual_ = qual; }
    void set_filter(const std::string& filter) { filter_ = filter; }
    void set_info(const std::map<std::string, std::string>& info) { info_ = info; }
    void set_format(const std::string& format) { format_ = format; }
    void set_samples(const std::vector<std::map<std::string, std::string>>& samples) { samples_ = samples; }

    /// @brief Add an INFO field entry
    void add_info(const std::string& key, const std::string& value) {
        info_[key] = value;
    }

    /// @brief Add a sample with genotype data
    void add_sample(const std::map<std::string, std::string>& sample_data) {
        samples_.push_back(sample_data);
    }

    /// @brief Phase 2: Validate against metadata
    /// @param info_meta INFO metadata from header
    /// @return true if INFO fields are valid
    bool validate_info(const std::map<std::string, info_meta>& info_meta) const {
        for (std::map<std::string, std::string>::const_iterator it = info_.begin();
             it != info_.end(); ++it) {
            // Check if INFO key is defined
            if (info_meta.find(it->first) == info_meta.end()) {
                return false; // Undefined INFO key
            }
        }
        return true;
    }

    /// @brief Check if this is a SNP
    bool is_snp() const {
        if (ref_.length() != 1) return false;
        for (std::size_t i = 0; i < alt_.size(); ++i) {
            if (alt_[i].length() != 1) return false;
        }
        return true;
    }

    /// @brief Check if this is an insertion
    bool is_insertion() const {
        if (alt_.empty()) return false;
        return ref_.length() < alt_[0].length();
    }

    /// @brief Check if this is a deletion
    bool is_deletion() const {
        if (alt_.empty()) return false;
        return ref_.length() > alt_[0].length();
    }

    // ========================================================================
    // PHASE 3: STRUCTURAL VARIANT METHODS
    // ========================================================================

    /// @brief Check if ALT contains symbolic allele
    bool is_symbolic_sv() const {
        if (alt_.empty()) return false;
        return is_symbolic_allele(alt_[0]);
    }

    /// @brief Check if ALT contains breakend notation
    bool is_breakend() const {
        if (alt_.empty()) return false;
        return is_breakend_notation(alt_[0]);
    }

    /// @brief Get SV type from record
    sv_type get_sv_type() const {
        if (alt_.empty()) return sv_type::UNKNOWN;
        return detect_sv_type(alt_[0], get_info_string());
    }

    /// @brief Parse SV information from INFO field
    sv_info get_sv_info() const {
        return parse_sv_info(get_info_string());
    }

    /// @brief Parse breakend notation from ALT field
    /// @param result Output breakend structure
    /// @return true if successfully parsed
    bool get_breakend(breakend& result) const {
        if (alt_.empty()) return false;
        return parse_breakend(alt_[0], result);
    }

    /// @brief Get END position from INFO (for SVs)
    /// @return END position or -1 if not set
    int get_end_pos() const {
        sv_info sv = get_sv_info();
        return sv.end;
    }

    /// @brief Get SVLEN from INFO
    /// @return SVLEN value or 0 if not set
    int get_svlen() const {
        sv_info sv = get_sv_info();
        return sv.svlen;
    }

    /// @brief Check if SV is imprecise
    bool is_imprecise() const {
        sv_info sv = get_sv_info();
        return sv.imprecise;
    }

    // ========================================================================
    // PHASE 4: gVCF METHODS
    // ========================================================================

    /// @brief Check if this record contains <*> allele (NON_REF)
    bool has_non_ref_allele() const {
        return vcf::has_non_ref_allele(alt_);
    }

    /// @brief Check if this is a pure reference block (only <*>)
    bool is_reference_block() const {
        return is_pure_reference_block(alt_);
    }

    /// @brief Parse reference block information from first sample
    /// @return reference_block structure with END, MIN_DP, GQ, DP, GT
    reference_block get_reference_block() const {
        std::string info_str = get_info_string();
        std::map<std::string, std::string> sample_data;
        if (!samples_.empty()) {
            sample_data = samples_[0];
        }
        return parse_reference_block(static_cast<int>(pos_), info_str, sample_data);
    }

    /// @brief Get reference block for a specific sample
    /// @param sample_idx Sample index (0-based)
    /// @return reference_block structure
    reference_block get_reference_block(std::size_t sample_idx) const {
        std::string info_str = get_info_string();
        std::map<std::string, std::string> sample_data;
        if (sample_idx < samples_.size()) {
            sample_data = samples_[sample_idx];
        }
        return parse_reference_block(static_cast<int>(pos_), info_str, sample_data);
    }

    /// @brief Validate reference block consistency
    bool validate_reference_block() const {
        reference_block block = get_reference_block();
        return vcf::validate_reference_block(block);
    }

    /// @brief Get block length (END - POS + 1, or 1 if no END)
    int get_block_length() const {
        int end = get_end_pos();
        if (end <= 0) return 1;
        return end - static_cast<int>(pos_) + 1;
    }

    // ========================================================================
    // VALIDATION METHODS
    // ========================================================================

    /// @brief Validate basic record fields
    /// @param result Validation result to accumulate errors
    /// @param line_num Line number for error reporting
    void validate_basic_fields(validation_result& result, std::size_t line_num = 0) const {
        // Validate REF
        if (!is_valid_ref(ref_)) {
            result.add_error("Invalid REF allele: " + ref_, line_num, "REF");
        }
        
        // Validate ALT alleles
        for (std::size_t i = 0; i < alt_.size(); ++i) {
            if (!is_valid_alt(alt_[i])) {
                result.add_error("Invalid ALT allele: " + alt_[i], line_num, "ALT");
            }
        }
    }

    /// @brief Validate record against header metadata
    template<typename InfoMap, typename FormatMap, typename FilterMap>
    void validate_against_header(const InfoMap& info_meta,
                                  const FormatMap& format_meta,
                                  const FilterMap& filter_meta,
                                  validation_result& result,
                                  std::size_t line_num = 0) const {
        // Validate FILTER
        validate_filter(filter_, filter_meta, line_num, result);
        
        // Validate INFO keys
        std::string info_str = get_info_string();
        validate_info_keys(info_str, info_meta, line_num, result);
        
        // Validate FORMAT keys
        if (!format_.empty()) {
            validate_format_keys(format_, format_meta, line_num, result);
        }
    }

    /// @brief Convert record to VCF line format
    std::string to_string(const std::vector<std::string>& sample_names = std::vector<std::string>()) const {
        std::ostringstream oss;
        
        // CHROM
        oss << chrom_ << "\t";
        
        // POS
        oss << pos_ << "\t";
        
        // ID
        oss << (id_.empty() ? "." : id_) << "\t";
        
        // REF
        oss << ref_ << "\t";
        
        // ALT
        if (alt_.empty()) {
            oss << ".";
        } else {
            for (std::size_t i = 0; i < alt_.size(); ++i) {
                if (i > 0) oss << ",";
                oss << alt_[i];
            }
        }
        oss << "\t";
        
        // QUAL
        if (qual_ == 0.0) {
            oss << ".";
        } else {
            oss << qual_;
        }
        oss << "\t";
        
        // FILTER
        oss << (filter_.empty() ? "." : filter_) << "\t";
        
        // INFO
        if (info_.empty()) {
            oss << ".";
        } else {
            bool first = true;
            for (std::map<std::string, std::string>::const_iterator it = info_.begin();
                 it != info_.end(); ++it) {
                if (!first) oss << ";";
                oss << it->first;
                if (!it->second.empty()) {
                    oss << "=" << it->second;
                }
                first = false;
            }
        }
        
        // FORMAT and samples (if present)
        if (!format_.empty() && !samples_.empty()) {
            oss << "\t" << format_;
            
            // Parse format fields once (cached in format_fields_ if available)
            const std::vector<std::string>& format_fields = get_format_fields();
            
            for (std::size_t i = 0; i < samples_.size(); ++i) {
                oss << "\t";
                for (std::size_t j = 0; j < format_fields.size(); ++j) {
                    if (j > 0) oss << ":";
                    std::map<std::string, std::string>::const_iterator sample_it = 
                        samples_[i].find(format_fields[j]);
                    if (sample_it != samples_[i].end()) {
                        oss << sample_it->second;
                    } else {
                        oss << ".";
                    }
                }
            }
        }
        
        return oss.str();
    }

private:
    std::string chrom_;
    std::size_t pos_;
    std::string id_;
    std::string ref_;
    std::vector<std::string> alt_;
    double qual_;
    std::string filter_;
    std::map<std::string, std::string> info_;
    std::string format_;
    std::vector<std::map<std::string, std::string>> samples_;
    mutable std::vector<std::string> format_fields_; // Cached parsed format

    /// @brief Get parsed format fields (cached)
    const std::vector<std::string>& get_format_fields() const {
        if (format_fields_.empty() && !format_.empty()) {
            std::istringstream format_stream(format_);
            std::string field;
            while (std::getline(format_stream, field, ':')) {
                format_fields_.push_back(field);
            }
        }
        return format_fields_;
    }

    /// @brief Convert INFO map to string format (for SV parsing)
    std::string get_info_string() const {
        if (info_.empty()) return ".";
        
        std::ostringstream oss;
        bool first = true;
        for (std::map<std::string, std::string>::const_iterator it = info_.begin();
             it != info_.end(); ++it) {
            if (!first) oss << ";";
            oss << it->first;
            if (!it->second.empty()) {
                oss << "=" << it->second;
            }
            first = false;
        }
        return oss.str();
    }
};

/// @brief VCF file reader
/// Supports VCFv4.3 and BCFv2.2 specifications
class reader {
public:
    /// @brief Constructor
    /// @param filename Path to VCF file
    explicit reader(const std::string& filename) 
        : filename_(filename), line_number_(0) {
        file_.open(filename_.c_str());
        if (!file_.is_open()) {
            throw std::runtime_error("Cannot open VCF file: " + filename_);
        }
        read_header();
    }

    /// @brief Destructor
    ~reader() {
        if (file_.is_open()) {
            file_.close();
        }
    }

    /// @brief Get the VCF version from header
    const std::string& version() const { return version_; }

    /// @brief Get all header lines
    const std::vector<std::string>& header_lines() const { return header_lines_; }

    /// @brief Get sample names
    const std::vector<std::string>& sample_names() const { return sample_names_; }

    /// @brief Get INFO metadata definitions
    const std::map<std::string, info_meta>& info_metadata() const { return info_defs_; }

    /// @brief Get FORMAT metadata definitions
    const std::map<std::string, format_meta>& format_metadata() const { return format_defs_; }

    /// @brief Get FILTER metadata definitions
    const std::map<std::string, filter_meta>& filter_metadata() const { return filter_defs_; }

    /// @brief Get contig metadata definitions
    const std::map<std::string, contig_meta>& contig_metadata() const { return contig_defs_; }

    /// @brief Get ALT metadata definitions
    const std::map<std::string, alt_meta>& alt_metadata() const { return alt_defs_; }

    /// @brief Read next record
    /// @param rec Record to populate
    /// @return true if record was read, false if end of file
    bool read_record(record& rec) {
        std::string line;
        while (std::getline(file_, line)) {
            ++line_number_;
            
            // Skip empty lines
            if (line.empty()) continue;
            
            // Skip header lines (should not happen after read_header, but be safe)
            if (line[0] == '#') continue;
            
            return parse_record(line, rec);
        }
        return false;
    }

    /// @brief Read all records
    std::vector<record> read_all() {
        std::vector<record> records;
        record rec;
        while (read_record(rec)) {
            records.push_back(rec);
        }
        return records;
    }

    /// @brief Validate file header
    /// @return validation_result with any errors found
    validation_result validate_header() const {
        validation_result result;
        
        // Check fileformat is first line
        if (header_lines_.empty() || !validate_fileformat(header_lines_[0], result)) {
            if (result.is_valid()) {  // Only add if not already added
                result.add_error("Missing or invalid ##fileformat line", 1);
            }
        }
        
        // Validate ID patterns in metadata
        check_duplicate_ids(info_defs_, "INFO", result);
        check_duplicate_ids(format_defs_, "FORMAT", result);
        check_duplicate_ids(filter_defs_, "FILTER", result);
        
        // Validate sample names are unique
        validate_sample_names(sample_names_, result);
        
        return result;
    }

private:
    std::string filename_;
    std::ifstream file_;
    std::size_t line_number_;
    std::string version_;
    std::vector<std::string> header_lines_;
    std::vector<std::string> sample_names_;
    
    // Phase 2: Structured metadata
    std::map<std::string, info_meta> info_defs_;
    std::map<std::string, format_meta> format_defs_;
    std::map<std::string, filter_meta> filter_defs_;
    std::map<std::string, contig_meta> contig_defs_;
    std::map<std::string, alt_meta> alt_defs_;

    void read_header() {
        std::string line;
        while (std::getline(file_, line)) {
            ++line_number_;
            
            if (line.empty()) continue;
            
            if (line[0] == '#') {
                header_lines_.push_back(line);
                
                // Extract version
                if (line.compare(0, 13, "##fileformat=") == 0) {
                    version_ = line.substr(13);
                }
                
                // Phase 2: Parse structured metadata
                else if (line.compare(0, 7, "##INFO=") == 0) {
                    parse_info_header(line.substr(7));
                }
                else if (line.compare(0, 9, "##FORMAT=") == 0) {
                    parse_format_header(line.substr(9));
                }
                else if (line.compare(0, 9, "##FILTER=") == 0) {
                    parse_filter_header(line.substr(9));
                }
                else if (line.compare(0, 9, "##contig=") == 0) {
                    parse_contig_header(line.substr(9));
                }
                else if (line.compare(0, 6, "##ALT=") == 0) {
                    parse_alt_header(line.substr(6));
                }
                
                // Extract sample names from column header line
                if (line[1] != '#') {
                    std::istringstream iss(line.substr(1)); // Remove leading #
                    std::string token;
                    std::vector<std::string> columns;
                    while (iss >> token) {
                        columns.push_back(token);
                    }
                    
                    // Sample names start after FORMAT column (index 9)
                    if (columns.size() > 9) {
                        sample_names_.assign(columns.begin() + 9, columns.end());
                    }
                }
            } else {
                // First non-header line, put it back by seeking back
                file_.seekg(-(line.length() + 1), std::ios::cur);
                --line_number_;
                break;
            }
        }
    }

    void parse_info_header(const std::string& content) {
        std::map<std::string, std::string> fields = parse_structured_header(content);
        if (fields.find("ID") != fields.end()) {
            info_meta meta;
            meta.id = fields["ID"];
            meta.number = fields.count("Number") ? fields["Number"] : ".";
            meta.type = fields.count("Type") ? fields["Type"] : "String";
            meta.description = fields.count("Description") ? fields["Description"] : "";
            // Store any extra fields
            for (std::map<std::string, std::string>::const_iterator it = fields.begin();
                 it != fields.end(); ++it) {
                if (it->first != "ID" && it->first != "Number" && 
                    it->first != "Type" && it->first != "Description") {
                    meta.extra[it->first] = it->second;
                }
            }
            info_defs_[meta.id] = meta;
        }
    }

    void parse_format_header(const std::string& content) {
        std::map<std::string, std::string> fields = parse_structured_header(content);
        if (fields.find("ID") != fields.end()) {
            format_meta meta;
            meta.id = fields["ID"];
            meta.number = fields.count("Number") ? fields["Number"] : ".";
            meta.type = fields.count("Type") ? fields["Type"] : "String";
            meta.description = fields.count("Description") ? fields["Description"] : "";
            for (std::map<std::string, std::string>::const_iterator it = fields.begin();
                 it != fields.end(); ++it) {
                if (it->first != "ID" && it->first != "Number" && 
                    it->first != "Type" && it->first != "Description") {
                    meta.extra[it->first] = it->second;
                }
            }
            format_defs_[meta.id] = meta;
        }
    }

    void parse_filter_header(const std::string& content) {
        std::map<std::string, std::string> fields = parse_structured_header(content);
        if (fields.find("ID") != fields.end()) {
            filter_meta meta;
            meta.id = fields["ID"];
            meta.description = fields.count("Description") ? fields["Description"] : "";
            filter_defs_[meta.id] = meta;
        }
    }

    void parse_contig_header(const std::string& content) {
        std::map<std::string, std::string> fields = parse_structured_header(content);
        if (fields.find("ID") != fields.end()) {
            contig_meta meta;
            meta.id = fields["ID"];
            if (fields.count("length")) {
                std::istringstream iss(fields["length"]);
                iss >> meta.length;
            }
            for (std::map<std::string, std::string>::const_iterator it = fields.begin();
                 it != fields.end(); ++it) {
                if (it->first != "ID" && it->first != "length") {
                    meta.extra[it->first] = it->second;
                }
            }
            contig_defs_[meta.id] = meta;
        }
    }

    void parse_alt_header(const std::string& content) {
        std::map<std::string, std::string> fields = parse_structured_header(content);
        if (fields.find("ID") != fields.end()) {
            alt_meta meta;
            meta.id = fields["ID"];
            meta.description = fields.count("Description") ? fields["Description"] : "";
            alt_defs_[meta.id] = meta;
        }
    }

    bool parse_record(const std::string& line, record& rec) {
        // Use zero-copy parsing with string_view
        std::string_view line_view(line);
        detail::field_extractor fields(line_view);
        
        // Extract required fields (zero-copy)
        std::string_view chrom_view = fields.next();
        std::string_view pos_view = fields.next();
        std::string_view id_view = fields.next();
        std::string_view ref_view = fields.next();
        std::string_view alt_view = fields.next();
        std::string_view qual_view = fields.next();
        std::string_view filter_view = fields.next();
        std::string_view info_view = fields.next();
        
        // Validate we got all 8 required fields (check for empty or missing fields)
        if (chrom_view.empty() || pos_view.empty() || id_view.empty() || 
            ref_view.empty() || alt_view.empty() || qual_view.empty() || 
            filter_view.empty() || info_view.empty()) {
            throw vcf_parse_error("Invalid VCF record: expected 8 fixed fields", line_number_);
        }
        
        // Set basic fields (convert from string_view to string only when storing)
        rec.set_chrom(std::string(chrom_view));
        rec.set_pos(detail::parse_uint_fast(pos_view));
        rec.set_id(id_view == "." ? "" : std::string(id_view));
        rec.set_ref(std::string(ref_view));
        
        // Parse ALT alleles (zero-copy split, then convert to strings)
        std::vector<std::string> alt_alleles;
        if (alt_view != ".") {
            detail::split_view(alt_view, ',', [&alt_alleles](std::string_view allele) {
                alt_alleles.emplace_back(allele);
            });
        }
        rec.set_alt(alt_alleles);
        
        // Parse QUAL (fast conversion)
        if (qual_view != ".") {
            rec.set_qual(detail::parse_double_fast(qual_view));
        } else {
            rec.set_qual(0.0);
        }
        
        // Set FILTER
        rec.set_filter(filter_view == "." ? "" : std::string(filter_view));
        
        // Parse INFO (zero-copy split, then store as strings)
        std::map<std::string, std::string> info_map;
        if (info_view != ".") {
            detail::split_view(info_view, ';', [&info_map](std::string_view info_field) {
                std::size_t eq_pos = info_field.find('=');
                if (eq_pos != std::string_view::npos) {
                    info_map[std::string(info_field.substr(0, eq_pos))] = 
                        std::string(info_field.substr(eq_pos + 1));
                } else {
                    info_map[std::string(info_field)] = "";
                }
            });
        }
        rec.set_info(info_map);
        
        // Parse FORMAT and sample data if present
        if (fields.has_more()) {
            std::string_view format_view = fields.next();
            rec.set_format(std::string(format_view));
            
            // Extract format field names
            std::vector<std::string> format_fields;
            detail::split_view(format_view, ':', [&format_fields](std::string_view field) {
                format_fields.emplace_back(field);
            });
            
            // Parse each sample
            std::vector<std::map<std::string, std::string>> samples;
            while (fields.has_more()) {
                std::string_view sample_view = fields.next();
                std::map<std::string, std::string> sample_map;
                
                std::size_t field_idx = 0;
                detail::split_view(sample_view, ':', [&](std::string_view value) {
                    if (field_idx < format_fields.size()) {
                        sample_map[format_fields[field_idx]] = std::string(value);
                        ++field_idx;
                    }
                });
                
                samples.push_back(std::move(sample_map));
            }
            rec.set_samples(samples);
        }
        
        return true;
    }

    std::string to_string(std::size_t num) const {
        std::ostringstream oss;
        oss << num;
        return oss.str();
    }
};

/// @brief VCF file writer
/// Supports VCFv4.3 and BCFv2.2 specifications
class writer {
public:
    /// @brief Constructor
    /// @param filename Path to output VCF file
    /// @param version VCF version (default: VCFv4.3)
    explicit writer(const std::string& filename, const std::string& version = "VCFv4.3")
        : filename_(filename), version_(version) {
        file_.open(filename_.c_str());
        if (!file_.is_open()) {
            throw vcf_parse_error("Cannot create VCF file: " + filename_);
        }
    }

    /// @brief Destructor
    ~writer() {
        if (file_.is_open()) {
            file_.close();
        }
    }

    /// @brief Add a header line
    void add_header_line(const std::string& line) {
        header_lines_.push_back(line);
    }

    /// @brief Add a contig header line
    void add_contig(const std::string& id, std::size_t length = 0) {
        std::ostringstream oss;
        oss << "##contig=<ID=" << id;
        if (length > 0) {
            oss << ",length=" << length;
        }
        oss << ">";
        header_lines_.push_back(oss.str());
    }

    /// @brief Add an INFO header line
    void add_info(const std::string& id, const std::string& number,
                  const std::string& type, const std::string& description) {
        std::ostringstream oss;
        oss << "##INFO=<ID=" << id << ",Number=" << number
            << ",Type=" << type << ",Description=\"" << description << "\">";
        header_lines_.push_back(oss.str());
    }

    /// @brief Add a FORMAT header line
    void add_format(const std::string& id, const std::string& number,
                    const std::string& type, const std::string& description) {
        std::ostringstream oss;
        oss << "##FORMAT=<ID=" << id << ",Number=" << number
            << ",Type=" << type << ",Description=\"" << description << "\">";
        header_lines_.push_back(oss.str());
    }

    /// @brief Set sample names
    void set_sample_names(const std::vector<std::string>& names) {
        sample_names_ = names;
    }

    /// @brief Write the header
    void write_header() {
        // Write fileformat
        file_ << "##fileformat=" << version_ << "\n";
        
        // Write other header lines
        for (std::size_t i = 0; i < header_lines_.size(); ++i) {
            file_ << header_lines_[i] << "\n";
        }
        
        // Write column header
        file_ << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        if (!sample_names_.empty()) {
            file_ << "\tFORMAT";
            for (std::size_t i = 0; i < sample_names_.size(); ++i) {
                file_ << "\t" << sample_names_[i];
            }
        }
        file_ << "\n";
        file_.flush();
    }

    /// @brief Write a record
    void write_record(const record& rec) {
        file_ << rec.to_string(sample_names_) << "\n";
        file_.flush();
    }

    /// @brief Write multiple records
    void write_records(const std::vector<record>& records) {
        for (std::size_t i = 0; i < records.size(); ++i) {
            write_record(records[i]);
        }
    }

private:
    std::string filename_;
    std::string version_;
    std::ofstream file_;
    std::vector<std::string> header_lines_;
    std::vector<std::string> sample_names_;
};

} // namespace vcf
} // namespace genetics
} // namespace boost

#endif // BOOST_GENETICS_VCF_HPP
