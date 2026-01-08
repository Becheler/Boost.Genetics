#ifndef BOOST_GENETICS_BCF_HPP
#define BOOST_GENETICS_BCF_HPP

#include <string>
#include <vector>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <fstream>

// Conditional Boost.IOStreams support
#ifdef BOOST_GENETICS_HAS_IOSTREAMS
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#endif

namespace boost {
namespace genetics {
namespace bcf {

// ============================================================================
// BCF BINARY FORMAT SUPPORT (VCF v4.3 / BCF v2.2)
// ============================================================================

/// BCF file magic bytes
static const char BCF_MAGIC[3] = {'B', 'C', 'F'};

/// BCF version constants
static const uint8_t BCF_MAJOR_VERSION = 2;
static const uint8_t BCF_MINOR_VERSION = 2;

/// BCF type codes (lower 4 bits of type descriptor)
enum class bcf_type : uint8_t {
    MISSING = 0x0,  ///< Missing value
    INT8    = 0x1,  ///< 8-bit integer
    INT16   = 0x2,  ///< 16-bit integer
    INT32   = 0x3,  ///< 32-bit integer
    FLOAT   = 0x5,  ///< 32-bit float
    CHAR    = 0x7   ///< Character/string
};

/// Special integer values for missing/end-of-vector
struct bcf_int_special {
    static const int8_t  INT8_MISSING  = static_cast<int8_t>(0x80);
    static const int8_t  INT8_EOV      = static_cast<int8_t>(0x81);
    static const int16_t INT16_MISSING = static_cast<int16_t>(0x8000);
    static const int16_t INT16_EOV     = static_cast<int16_t>(0x8001);
    static const int32_t INT32_MISSING = static_cast<int32_t>(0x80000000);
    static const int32_t INT32_EOV     = static_cast<int32_t>(0x80000001);
};

/// Special float values
struct bcf_float_special {
    static const uint32_t MISSING = 0x7F800001;
    static const uint32_t EOV     = 0x7F800002;
};

// ============================================================================
// TYPE DESCRIPTOR ENCODING
// ============================================================================

/// Type descriptor byte structure:
/// Bits 5-8 (upper 4 bits): Size (or 15 for overflow)
/// Bits 1-4 (lower 4 bits): Type code
inline uint8_t make_type_descriptor(bcf_type type, uint8_t size) {
    if (size >= 15) {
        return (15 << 4) | static_cast<uint8_t>(type);
    }
    return (size << 4) | static_cast<uint8_t>(type);
}

/// Extract type from descriptor
inline bcf_type get_type(uint8_t descriptor) {
    return static_cast<bcf_type>(descriptor & 0x0F);
}

/// Extract size from descriptor (returns 15 if overflow)
inline uint8_t get_size(uint8_t descriptor) {
    return (descriptor >> 4) & 0x0F;
}

// ============================================================================
// GENOTYPE ENCODING
// ============================================================================

/// Encode VCF genotype allele to BCF format: (allele + 1) << 1 | phased
/// @param allele Allele index (-1 for missing)
/// @param phased true if phased (|), false if unphased (/)
inline uint8_t encode_gt_allele(int allele, bool phased) {
    if (allele < 0) {
        return 0;  // Missing allele
    }
    return static_cast<uint8_t>(((allele + 1) << 1) | (phased ? 1 : 0));
}

/// Decode BCF genotype allele to VCF format
/// @param encoded Encoded allele byte
/// @param phased Output: true if phased
/// @return Allele index (-1 for missing)
inline int decode_gt_allele(uint8_t encoded, bool& phased) {
    if (encoded == 0) {
        phased = false;
        return -1;  // Missing
    }
    if (encoded == static_cast<uint8_t>(bcf_int_special::INT8_EOV)) {
        phased = false;
        return -2;  // End of vector
    }
    phased = (encoded & 1) != 0;
    return (encoded >> 1) - 1;
}

/// Encode VCF genotype string (e.g., "0/1", "1|0") to BCF bytes
inline std::vector<uint8_t> encode_genotype(const std::string& gt) {
    std::vector<uint8_t> result;
    
    if (gt.empty() || gt == ".") {
        result.push_back(0);
        return result;
    }
    
    std::string allele_str;
    bool is_phased = false;
    
    for (std::size_t i = 0; i < gt.size(); ++i) {
        char c = gt[i];
        
        if (c == '/' || c == '|') {
            // Process accumulated allele
            int allele = (allele_str == ".") ? -1 : std::atoi(allele_str.c_str());
            result.push_back(encode_gt_allele(allele, is_phased));
            
            // Next allele will be phased if separator is |
            is_phased = (c == '|');
            allele_str.clear();
        } else {
            allele_str += c;
        }
    }
    
    // Process final allele
    if (!allele_str.empty()) {
        int allele = (allele_str == ".") ? -1 : std::atoi(allele_str.c_str());
        result.push_back(encode_gt_allele(allele, is_phased));
    }
    
    return result;
}

/// Decode BCF genotype bytes to VCF string
inline std::string decode_genotype(const std::vector<uint8_t>& data) {
    if (data.empty()) return ".";
    
    std::string result;
    bool prev_phased = false;
    
    for (std::size_t i = 0; i < data.size(); ++i) {
        if (data[i] == static_cast<uint8_t>(bcf_int_special::INT8_EOV)) {
            break;  // End of vector
        }
        
        bool phased;
        int allele = decode_gt_allele(data[i], phased);
        
        if (i > 0) {
            // Use the phasing bit from the current (second) allele
            result += phased ? '|' : '/';
        }
        
        if (allele < 0) {
            result += '.';
        } else {
            result += std::to_string(allele);
        }
        
        prev_phased = phased;
    }
    
    return result.empty() ? "." : result;
}

// ============================================================================
// BCF HEADER
// ============================================================================

/// BCF file header structure
struct bcf_header {
    uint8_t major_version;   ///< Major version (2)
    uint8_t minor_version;   ///< Minor version (2)
    uint32_t l_text;         ///< Length of VCF header text
    std::string text;        ///< VCF header text (NUL-terminated)
    
    bcf_header() : major_version(BCF_MAJOR_VERSION), 
                   minor_version(BCF_MINOR_VERSION), 
                   l_text(0) {}
};

/// Read BCF header from binary stream
inline bool read_bcf_header(std::istream& in, bcf_header& header) {
    // Read magic bytes
    char magic[3];
    in.read(magic, 3);
    if (!in || std::memcmp(magic, BCF_MAGIC, 3) != 0) {
        return false;  // Not a BCF file
    }
    
    // Read version
    in.read(reinterpret_cast<char*>(&header.major_version), 1);
    in.read(reinterpret_cast<char*>(&header.minor_version), 1);
    
    if (!in) return false;
    
    // Validate version
    if (header.major_version != BCF_MAJOR_VERSION) {
        return false;  // Unsupported version
    }
    
    // Read text length (little-endian)
    in.read(reinterpret_cast<char*>(&header.l_text), sizeof(uint32_t));
    if (!in) return false;
    
    // Read text
    if (header.l_text > 0) {
        std::vector<char> buffer(header.l_text);
        in.read(buffer.data(), header.l_text);
        if (!in) return false;
        
        header.text.assign(buffer.begin(), buffer.end());
        
        // Remove trailing NUL if present
        if (!header.text.empty() && header.text.back() == '\0') {
            header.text.pop_back();
        }
    }
    
    return true;
}

/// Write BCF header to binary stream
inline bool write_bcf_header(std::ostream& out, const bcf_header& header) {
    // Write magic bytes
    out.write(BCF_MAGIC, 3);
    
    // Write version
    out.write(reinterpret_cast<const char*>(&header.major_version), 1);
    out.write(reinterpret_cast<const char*>(&header.minor_version), 1);
    
    // Write text length
    uint32_t l_text = static_cast<uint32_t>(header.text.size() + 1);  // +1 for NUL
    out.write(reinterpret_cast<const char*>(&l_text), sizeof(uint32_t));
    
    // Write text with NUL terminator
    out.write(header.text.c_str(), header.text.size());
    out.put('\0');
    
    return out.good();
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

/// Check if value is BCF missing marker
inline bool is_bcf_missing_int8(int8_t val) {
    return val == bcf_int_special::INT8_MISSING;
}

inline bool is_bcf_missing_int16(int16_t val) {
    return val == bcf_int_special::INT16_MISSING;
}

inline bool is_bcf_missing_int32(int32_t val) {
    return val == bcf_int_special::INT32_MISSING;
}

/// Check if value is end-of-vector marker
inline bool is_bcf_eov_int8(int8_t val) {
    return val == bcf_int_special::INT8_EOV;
}

inline bool is_bcf_eov_int16(int16_t val) {
    return val == bcf_int_special::INT16_EOV;
}

inline bool is_bcf_eov_int32(int32_t val) {
    return val == bcf_int_special::INT32_EOV;
}

#ifdef BOOST_GENETICS_HAS_IOSTREAMS
// ============================================================================
// BGZF COMPRESSED BCF READER (requires Boost.IOStreams)
// ============================================================================

/// BCF reader with BGZF decompression
class bcf_reader {
public:
    /// Constructor
    /// @param filename Path to BCF file (BGZF compressed)
    explicit bcf_reader(const std::string& filename) {
        // Create filtering stream with gzip decompressor
        in_.push(boost::iostreams::gzip_decompressor());
        in_.push(boost::iostreams::file_source(filename, std::ios::binary));
        
        // Read header
        if (!read_bcf_header(in_, header_)) {
            throw std::runtime_error("Failed to read BCF header from: " + filename);
        }
    }
    
    /// Get BCF header
    const bcf_header& header() const { return header_; }
    
    /// Get VCF header text
    const std::string& vcf_header_text() const { return header_.text; }
    
    /// Check if more records available
    bool has_more() const { return in_.good(); }
    
private:
    boost::iostreams::filtering_istream in_;
    bcf_header header_;
};

/// BCF writer with BGZF compression
class bcf_writer {
public:
    /// Constructor
    /// @param filename Path to output BCF file
    /// @param vcf_header VCF header text
    explicit bcf_writer(const std::string& filename, const std::string& vcf_header) {
        // Create filtering stream with gzip compressor
        out_.push(boost::iostreams::gzip_compressor());
        out_.push(boost::iostreams::file_sink(filename, std::ios::binary));
        
        // Write header
        bcf_header header;
        header.text = vcf_header;
        if (!write_bcf_header(out_, header)) {
            throw std::runtime_error("Failed to write BCF header to: " + filename);
        }
    }
    
    /// Destructor - flushes and closes stream
    ~bcf_writer() {
        out_.reset();
    }
    
private:
    boost::iostreams::filtering_ostream out_;
};

#endif // BOOST_GENETICS_HAS_IOSTREAMS

} // namespace bcf
} // namespace genetics
} // namespace boost

#endif // BOOST_GENETICS_BCF_HPP
