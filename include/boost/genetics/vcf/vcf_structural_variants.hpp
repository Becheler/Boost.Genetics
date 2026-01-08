#ifndef BOOST_GENETICS_VCF_STRUCTURAL_VARIANTS_HPP
#define BOOST_GENETICS_VCF_STRUCTURAL_VARIANTS_HPP

#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <stdexcept>

namespace boost {
namespace genetics {
namespace vcf {

// ============================================================================
// PHASE 3: STRUCTURAL VARIANTS
// ============================================================================

/// Structural variant types as defined in VCF v4.3
enum class sv_type {
    DEL,      ///< Deletion
    INS,      ///< Insertion
    DUP,      ///< Duplication
    INV,      ///< Inversion
    CNV,      ///< Copy number variation
    BND,      ///< Breakend
    UNKNOWN   ///< Unknown or not an SV
};

/// Convert SV type enum to string
inline std::string sv_type_to_string(sv_type type) {
    switch (type) {
        case sv_type::DEL: return "DEL";
        case sv_type::INS: return "INS";
        case sv_type::DUP: return "DUP";
        case sv_type::INV: return "INV";
        case sv_type::CNV: return "CNV";
        case sv_type::BND: return "BND";
        case sv_type::UNKNOWN: return "UNKNOWN";
        default: return "UNKNOWN";
    }
}

/// Convert string to SV type enum
inline sv_type string_to_sv_type(const std::string& str) {
    if (str == "DEL") return sv_type::DEL;
    if (str == "INS") return sv_type::INS;
    if (str == "DUP") return sv_type::DUP;
    if (str == "INV") return sv_type::INV;
    if (str == "CNV") return sv_type::CNV;
    if (str == "BND") return sv_type::BND;
    return sv_type::UNKNOWN;
}

// ============================================================================
// BREAKEND NOTATION
// ============================================================================

/// Breakend orientation
enum class breakend_orientation {
    FORWARD,  ///< Forward strand
    REVERSE   ///< Reverse complement strand
};

/// Breakend position (left or right of reference base)
enum class breakend_position {
    BEFORE,   ///< Before the reference base(s)
    AFTER     ///< After the reference base(s)
};

/// Parsed breakend notation
struct breakend {
    std::string ref_bases;              ///< Reference bases being replaced
    std::string novel_sequence;          ///< Novel inserted sequence (if any)
    std::string mate_chr;                ///< Chromosome of mate breakend
    int mate_pos;                        ///< Position of mate breakend
    breakend_orientation orientation;    ///< Strand orientation
    breakend_position position;          ///< Whether mate extends before/after
    
    breakend() : mate_pos(0), orientation(breakend_orientation::FORWARD), 
                 position(breakend_position::AFTER) {}
};

/// Parse breakend notation from ALT field
/// Formats: t[p[, t]p], ]p]t, [p[t
/// Returns true if successfully parsed as a breakend
inline bool parse_breakend(std::string_view alt, breakend& result) {
    // Try to find bracket patterns
    std::size_t open_bracket = alt.find('[');
    std::size_t close_bracket = alt.find(']');
    
    if (open_bracket == std::string_view::npos && close_bracket == std::string_view::npos) {
        return false;  // Not a breakend
    }
    
    // Pattern: t[p[ or [p[t
    if (open_bracket != std::string_view::npos) {
        std::size_t second_bracket = alt.find('[', open_bracket + 1);
        if (second_bracket == std::string_view::npos) return false;
        
        // Extract position between brackets
        std::string_view pos_str = alt.substr(open_bracket + 1, second_bracket - open_bracket - 1);
        
        // Split chr:pos
        std::size_t colon = pos_str.find(':');
        if (colon == std::string_view::npos) return false;
        
        result.mate_chr = std::string(pos_str.substr(0, colon));
        std::string pos_part(pos_str.substr(colon + 1));
        result.mate_pos = std::atoi(pos_part.c_str());
        
        // Determine orientation and position
        if (open_bracket == 0) {
            // [p[t - reverse complement extending right, joined before t
            result.novel_sequence = std::string(alt.substr(second_bracket + 1));
            result.orientation = breakend_orientation::REVERSE;
            result.position = breakend_position::BEFORE;
        } else {
            // t[p[ - piece extending right of p joined after t
            result.novel_sequence = std::string(alt.substr(0, open_bracket));
            result.orientation = breakend_orientation::FORWARD;
            result.position = breakend_position::AFTER;
        }
        
        return true;
    }
    
    // Pattern: t]p] or ]p]t
    if (close_bracket != std::string_view::npos) {
        std::size_t second_bracket = alt.find(']', close_bracket + 1);
        if (second_bracket == std::string_view::npos) return false;
        
        // Extract position between brackets
        std::string_view pos_str = alt.substr(close_bracket + 1, second_bracket - close_bracket - 1);
        
        // Split chr:pos
        std::size_t colon = pos_str.find(':');
        if (colon == std::string_view::npos) return false;
        
        result.mate_chr = std::string(pos_str.substr(0, colon));
        std::string pos_part(pos_str.substr(colon + 1));
        result.mate_pos = std::atoi(pos_part.c_str());
        
        // Determine orientation and position
        if (close_bracket == 0) {
            // ]p]t - piece extending left of p joined before t
            result.novel_sequence = std::string(alt.substr(second_bracket + 1));
            result.orientation = breakend_orientation::FORWARD;
            result.position = breakend_position::BEFORE;
        } else {
            // t]p] - reverse complement extending left of p joined after t
            result.novel_sequence = std::string(alt.substr(0, close_bracket));
            result.orientation = breakend_orientation::REVERSE;
            result.position = breakend_position::AFTER;
        }
        
        return true;
    }
    
    return false;
}

// ============================================================================
// STRUCTURAL VARIANT INFO
// ============================================================================

/// Structural variant information extracted from INFO field
struct sv_info {
    sv_type type;               ///< SV type
    bool imprecise;             ///< IMPRECISE flag
    bool novel;                 ///< NOVEL flag
    int end;                    ///< END position (-1 if not set)
    int svlen;                  ///< SVLEN value (0 if not set)
    std::pair<int, int> cipos;  ///< CIPOS confidence interval
    std::pair<int, int> ciend;  ///< CIEND confidence interval
    std::string event;          ///< EVENT ID
    std::vector<std::string> mateid;  ///< MATEID values
    
    sv_info() : type(sv_type::UNKNOWN), imprecise(false), novel(false), 
                end(-1), svlen(0), cipos(0, 0), ciend(0, 0) {}
};

/// Parse INFO field to extract SV-specific information
inline sv_info parse_sv_info(std::string_view info_str) {
    sv_info result;
    
    if (info_str.empty() || info_str == ".") {
        return result;
    }
    
    // Split by semicolon (convert to string for istringstream)
    std::string info_s(info_str);
    std::istringstream iss(info_s);
    std::string token;
    
    while (std::getline(iss, token, ';')) {
        if (token.empty()) continue;
        
        std::size_t eq = token.find('=');
        
        if (eq == std::string::npos) {
            // Flag field
            if (token == "IMPRECISE") result.imprecise = true;
            else if (token == "NOVEL") result.novel = true;
        } else {
            // Key=Value field
            std::string key = token.substr(0, eq);
            std::string value = token.substr(eq + 1);
            
            if (key == "SVTYPE") {
                result.type = string_to_sv_type(value);
            } else if (key == "END") {
                result.end = std::atoi(value.c_str());
            } else if (key == "SVLEN") {
                result.svlen = std::atoi(value.c_str());
            } else if (key == "EVENT") {
                result.event = value;
            } else if (key == "MATEID") {
                // Can be comma-separated list
                std::istringstream mate_iss(value);
                std::string mate;
                while (std::getline(mate_iss, mate, ',')) {
                    result.mateid.push_back(mate);
                }
            } else if (key == "CIPOS") {
                // Parse two integers
                std::istringstream ci_iss(value);
                std::string first, second;
                if (std::getline(ci_iss, first, ',') && std::getline(ci_iss, second, ',')) {
                    result.cipos = std::make_pair(std::atoi(first.c_str()), std::atoi(second.c_str()));
                }
            } else if (key == "CIEND") {
                // Parse two integers
                std::istringstream ci_iss(value);
                std::string first, second;
                if (std::getline(ci_iss, first, ',') && std::getline(ci_iss, second, ',')) {
                    result.ciend = std::make_pair(std::atoi(first.c_str()), std::atoi(second.c_str()));
                }
            }
        }
    }
    
    return result;
}

// ============================================================================
// SYMBOLIC ALLELE DETECTION
// ============================================================================

/// Check if ALT field is a symbolic allele
inline bool is_symbolic_allele(std::string_view alt) {
    return !alt.empty() && alt[0] == '<' && alt[alt.size() - 1] == '>';
}

/// Extract symbolic allele ID (e.g., "<DEL>" -> "DEL")
inline std::string get_symbolic_id(std::string_view alt) {
    if (!is_symbolic_allele(alt)) return "";
    return std::string(alt.substr(1, alt.size() - 2));
}

/// Check if ALT field is a breakend notation
inline bool is_breakend_notation(std::string_view alt) {
    return alt.find('[') != std::string_view::npos || alt.find(']') != std::string_view::npos;
}

// ============================================================================
// SV TYPE DETECTION
// ============================================================================

/// Determine SV type from ALT field and INFO
inline sv_type detect_sv_type(std::string_view alt, std::string_view info_str) {
    // Check symbolic alleles first
    if (is_symbolic_allele(alt)) {
        std::string id = get_symbolic_id(alt);
        return string_to_sv_type(id);
    }
    
    // Check breakend notation
    if (is_breakend_notation(alt)) {
        return sv_type::BND;
    }
    
    // Check INFO field for SVTYPE
    sv_info info = parse_sv_info(info_str);
    return info.type;
}

} // namespace vcf
} // namespace genetics
} // namespace boost

#endif // BOOST_GENETICS_VCF_STRUCTURAL_VARIANTS_HPP
