#ifndef BOOST_GENETICS_VCF_GVCF_HPP
#define BOOST_GENETICS_VCF_GVCF_HPP

#include <string>
#include <map>
#include <vector>
#include <sstream>

namespace boost {
namespace genetics {
namespace vcf {

// ============================================================================
// PHASE 4: gVCF SUPPORT
// ============================================================================

/// Reference block information for gVCF format
/// gVCF uses <*> symbolic allele to indicate "any other allele"
/// and END INFO field to represent blocks of reference positions
struct reference_block {
    int start_pos;          ///< Start position (POS field)
    int end_pos;            ///< End position (END INFO field, -1 if not set)
    int min_depth;          ///< Minimum depth in block (MIN_DP FORMAT field, -1 if not set)
    int genotype_quality;   ///< Genotype quality (GQ FORMAT field, -1 if not set)
    int depth;              ///< Average depth (DP FORMAT field, -1 if not set)
    std::string genotype;   ///< Genotype (GT FORMAT field)
    
    reference_block() : start_pos(0), end_pos(-1), min_depth(-1), 
                        genotype_quality(-1), depth(-1) {}
};

/// Check if ALT field contains the NON_REF symbolic allele <*>
inline bool is_non_ref_allele(std::string_view alt) {
    return alt == "<*>" || alt == "<NON_REF>";
}

/// Check if ALT vector contains only <*> (pure reference block)
inline bool is_pure_reference_block(const std::vector<std::string>& alt) {
    if (alt.empty()) return false;
    if (alt.size() > 1) return false;
    return is_non_ref_allele(alt[0]);
}

/// Check if ALT vector contains <*> along with other alleles (mixed variant)
inline bool has_non_ref_allele(const std::vector<std::string>& alt) {
    for (std::size_t i = 0; i < alt.size(); ++i) {
        if (is_non_ref_allele(alt[i])) return true;
    }
    return false;
}

/// Parse reference block information from a record
/// Extracts END from INFO and MIN_DP, GQ, DP, GT from first sample
inline reference_block parse_reference_block(
    int pos,
    const std::string& info_str,
    const std::map<std::string, std::string>& sample_data)
{
    reference_block block;
    block.start_pos = pos;
    
    // Parse END from INFO field
    if (!info_str.empty() && info_str != ".") {
        std::istringstream iss(info_str);
        std::string token;
        while (std::getline(iss, token, ';')) {
            if (token.empty()) continue;
            std::size_t eq = token.find('=');
            if (eq != std::string::npos) {
                std::string key = token.substr(0, eq);
                std::string value = token.substr(eq + 1);
                if (key == "END") {
                    block.end_pos = std::atoi(value.c_str());
                }
            }
        }
    }
    
    // Parse FORMAT fields from sample data
    std::map<std::string, std::string>::const_iterator it;
    
    it = sample_data.find("MIN_DP");
    if (it != sample_data.end() && it->second != ".") {
        block.min_depth = std::atoi(it->second.c_str());
    }
    
    it = sample_data.find("GQ");
    if (it != sample_data.end() && it->second != ".") {
        block.genotype_quality = std::atoi(it->second.c_str());
    }
    
    it = sample_data.find("DP");
    if (it != sample_data.end() && it->second != ".") {
        block.depth = std::atoi(it->second.c_str());
    }
    
    it = sample_data.find("GT");
    if (it != sample_data.end()) {
        block.genotype = it->second;
    }
    
    return block;
}

/// Check if genotype is reference (0/0 or 0|0)
inline bool is_reference_genotype(const std::string& gt) {
    return gt == "0/0" || gt == "0|0" || gt == "0";
}

/// Check if genotype indicates no-call (./. or .|.)
inline bool is_no_call_genotype(const std::string& gt) {
    return gt == "./." || gt == ".|." || gt == ".";
}

/// Calculate block length from start and end positions
inline int get_block_length(const reference_block& block) {
    if (block.end_pos <= 0) return 1;  // Single position
    return block.end_pos - block.start_pos + 1;
}

/// Validate reference block consistency
/// Returns true if block is valid, false otherwise
inline bool validate_reference_block(const reference_block& block) {
    // Must have valid start position
    if (block.start_pos <= 0) return false;
    
    // If END is set, it must be >= start
    if (block.end_pos > 0 && block.end_pos < block.start_pos) return false;
    
    // Reference blocks should have reference genotype
    if (!block.genotype.empty() && !is_reference_genotype(block.genotype)) {
        // Allow no-call genotypes in some cases
        if (!is_no_call_genotype(block.genotype)) {
            return false;
        }
    }
    
    // MIN_DP should not exceed DP if both are set
    if (block.min_depth > 0 && block.depth > 0) {
        if (block.min_depth > block.depth) return false;
    }
    
    return true;
}

/// gVCF-specific FORMAT field definitions
struct gvcf_format_keys {
    static const char* MIN_DP;  ///< Minimum depth in reference block
    static const char* GQ;      ///< Genotype quality
    static const char* DP;      ///< Read depth
    static const char* PL;      ///< Phred-scaled genotype likelihoods
};

const char* gvcf_format_keys::MIN_DP = "MIN_DP";
const char* gvcf_format_keys::GQ = "GQ";
const char* gvcf_format_keys::DP = "DP";
const char* gvcf_format_keys::PL = "PL";

/// gVCF-specific INFO field definitions
struct gvcf_info_keys {
    static const char* END;  ///< End position of reference block
};

const char* gvcf_info_keys::END = "END";

} // namespace vcf
} // namespace genetics
} // namespace boost

#endif // BOOST_GENETICS_VCF_GVCF_HPP
