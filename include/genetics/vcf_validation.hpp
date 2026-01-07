#ifndef BOOST_GENETICS_VCF_VALIDATION_HPP
#define BOOST_GENETICS_VCF_VALIDATION_HPP

#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <cctype>
#include <cstdlib>

namespace boost {
namespace genetics {
namespace vcf {

// ============================================================================
// VALIDATION FUNCTIONS (Based on VCF v4.3 Specification)
// ============================================================================

/// Validation error information
struct validation_error {
    std::string message;    ///< Error description
    std::size_t line;       ///< Line number (0 if not applicable)
    std::string field;      ///< Field name where error occurred
    
    validation_error(const std::string& msg, std::size_t ln = 0, const std::string& fld = "")
        : message(msg), line(ln), field(fld) {}
};

/// Collection of validation errors
class validation_result {
public:
    validation_result() : valid_(true) {}
    
    /// Add an error
    void add_error(const std::string& msg, std::size_t line = 0, const std::string& field = "") {
        errors_.push_back(validation_error(msg, line, field));
        valid_ = false;
    }
    
    /// Check if validation passed
    bool is_valid() const { return valid_; }
    
    /// Get all errors
    const std::vector<validation_error>& errors() const { return errors_; }
    
    /// Get error count
    std::size_t error_count() const { return errors_.size(); }
    
private:
    bool valid_;
    std::vector<validation_error> errors_;
};

// ============================================================================
// ID PATTERN VALIDATION
// ============================================================================

/// Validate VCF ID pattern: ^([A-Za-z_][0-9A-Za-z_.]*)$
/// Exception: "1000G" is allowed despite starting with a digit
inline bool is_valid_vcf_id(const std::string& id) {
    if (id.empty()) return false;
    
    // Special case for "1000G"
    if (id == "1000G") return true;
    
    // First character must be letter or underscore
    if (!std::isalpha(id[0]) && id[0] != '_') return false;
    
    // Remaining characters: alphanumeric, underscore, or period
    for (std::size_t i = 1; i < id.size(); ++i) {
        char c = id[i];
        if (!std::isalnum(c) && c != '_' && c != '.') {
            return false;
        }
    }
    
    return true;
}

// ============================================================================
// POSITION VALIDATION
// ============================================================================

/// Validate POS field (must be positive integer, or 0 for telomere)
inline bool is_valid_pos(const std::string& pos_str) {
    if (pos_str.empty()) return false;
    
    // Check all digits
    for (std::size_t i = 0; i < pos_str.size(); ++i) {
        if (!std::isdigit(pos_str[i])) return false;
    }
    
    // Convert to integer and check >= 0
    long pos = std::atol(pos_str.c_str());
    return pos >= 0;
}

// ============================================================================
// ALLELE VALIDATION
// ============================================================================

/// Check if character is a valid nucleotide (ACGTN, case-insensitive)
inline bool is_nucleotide(char c) {
    char upper = std::toupper(c);
    return upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T' || upper == 'N';
}

/// Validate REF allele (cannot be empty, must be nucleotides)
inline bool is_valid_ref(const std::string& ref) {
    if (ref.empty()) return false;
    
    for (std::size_t i = 0; i < ref.size(); ++i) {
        if (!is_nucleotide(ref[i])) return false;
    }
    
    return true;
}

/// Validate ALT allele (nucleotides, *, ., or symbolic <ID>)
inline bool is_valid_alt(const std::string& alt) {
    if (alt.empty()) return false;
    
    // Dot notation (no variant)
    if (alt == ".") return true;
    
    // Star notation (overlapping deletion)
    if (alt == "*") return true;
    
    // Symbolic allele <ID>
    if (alt.size() >= 3 && alt[0] == '<' && alt[alt.size() - 1] == '>') {
        std::string id = alt.substr(1, alt.size() - 2);
        return is_valid_vcf_id(id);
    }
    
    // Breakend notation (contains [ or ])
    if (alt.find('[') != std::string::npos || alt.find(']') != std::string::npos) {
        return true;  // Complex validation handled separately
    }
    
    // Regular nucleotide sequence
    for (std::size_t i = 0; i < alt.size(); ++i) {
        if (!is_nucleotide(alt[i])) return false;
    }
    
    return true;
}

// ============================================================================
// QUALITY VALIDATION
// ============================================================================

/// Validate QUAL field (must be number or .)
inline bool is_valid_qual(const std::string& qual) {
    if (qual.empty()) return false;
    if (qual == ".") return true;
    
    // Simple float validation
    bool has_digit = false;
    bool has_dot = false;
    std::size_t start = 0;
    
    // Allow leading sign
    if (qual[0] == '-' || qual[0] == '+') start = 1;
    
    for (std::size_t i = start; i < qual.size(); ++i) {
        if (std::isdigit(qual[i])) {
            has_digit = true;
        } else if (qual[i] == '.') {
            if (has_dot) return false;  // Multiple dots
            has_dot = true;
        } else {
            return false;  // Invalid character
        }
    }
    
    return has_digit;
}

// ============================================================================
// HEADER VALIDATION
// ============================================================================

/// Validate that fileformat line is correct and first
inline bool validate_fileformat(const std::string& line, validation_result& result) {
    if (line.find("##fileformat=VCFv4.") != 0) {
        result.add_error("First line must be ##fileformat=VCFv4.x", 1);
        return false;
    }
    return true;
}

/// Check for duplicate IDs in metadata
template<typename MetaMap>
inline void check_duplicate_ids(const MetaMap& meta_map, const std::string& type, 
                                  validation_result& result) {
    // The map itself ensures uniqueness, but we can check if IDs are valid
    for (typename MetaMap::const_iterator it = meta_map.begin(); 
         it != meta_map.end(); ++it) {
        if (!is_valid_vcf_id(it->first)) {
            result.add_error(type + " ID '" + it->first + "' does not match required pattern", 
                            0, it->first);
        }
    }
}

/// Validate sample names are unique
inline void validate_sample_names(const std::vector<std::string>& samples, 
                                   validation_result& result) {
    std::set<std::string> seen;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (!seen.insert(samples[i]).second) {
            result.add_error("Duplicate sample name: " + samples[i], 0, "samples");
        }
    }
}

// ============================================================================
// RECORD VALIDATION
// ============================================================================

/// Validate that a record has correct number of fields
inline bool validate_field_count(const std::vector<std::string>& fields, 
                                  std::size_t expected_samples,
                                  std::size_t line_num,
                                  validation_result& result) {
    // Must have 8 fixed fields + optional FORMAT + samples
    if (fields.size() < 8) {
        result.add_error("Record has fewer than 8 required fields", line_num);
        return false;
    }
    
    if (expected_samples > 0) {
        // With samples: 8 fixed + 1 FORMAT + N samples = 9 + N
        std::size_t expected = 9 + expected_samples;
        if (fields.size() != expected) {
            std::ostringstream oss;
            oss << "Record has " << fields.size() << " fields, expected " << expected;
            result.add_error(oss.str(), line_num);
            return false;
        }
    } else {
        // Without samples: exactly 8 fields
        if (fields.size() != 8) {
            std::ostringstream oss;
            oss << "Record has " << fields.size() << " fields, expected 8";
            result.add_error(oss.str(), line_num);
            return false;
        }
    }
    
    return true;
}

/// Validate FILTER field against header definitions
template<typename FilterMap>
inline void validate_filter(const std::string& filter, 
                            const FilterMap& filter_defs,
                            std::size_t line_num,
                            validation_result& result) {
    if (filter.empty() || filter == "." || filter == "PASS") {
        return;  // Valid
    }
    
    // Can be semicolon-separated list
    std::istringstream iss(filter);
    std::string flt;
    while (std::getline(iss, flt, ';')) {
        if (flt != "PASS" && filter_defs.find(flt) == filter_defs.end()) {
            result.add_error("Undefined FILTER: " + flt, line_num, "FILTER");
        }
    }
}

/// Validate INFO keys against header definitions
template<typename InfoMap>
inline void validate_info_keys(const std::string& info_str,
                               const InfoMap& info_defs,
                               std::size_t line_num,
                               validation_result& result) {
    if (info_str.empty() || info_str == ".") return;
    
    std::istringstream iss(info_str);
    std::string token;
    while (std::getline(iss, token, ';')) {
        if (token.empty()) continue;
        
        std::size_t eq = token.find('=');
        std::string key = (eq == std::string::npos) ? token : token.substr(0, eq);
        
        if (info_defs.find(key) == info_defs.end()) {
            result.add_error("Undefined INFO key: " + key, line_num, "INFO");
        }
    }
}

/// Validate FORMAT keys against header definitions
template<typename FormatMap>
inline void validate_format_keys(const std::string& format_str,
                                 const FormatMap& format_defs,
                                 std::size_t line_num,
                                 validation_result& result) {
    if (format_str.empty()) return;
    
    std::istringstream iss(format_str);
    std::string key;
    while (std::getline(iss, key, ':')) {
        if (format_defs.find(key) == format_defs.end()) {
            result.add_error("Undefined FORMAT key: " + key, line_num, "FORMAT");
        }
    }
}

} // namespace vcf
} // namespace genetics
} // namespace boost

#endif // BOOST_GENETICS_VCF_VALIDATION_HPP
