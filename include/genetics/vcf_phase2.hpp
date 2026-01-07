//  Copyright (c) 2026 Boost.Genetics
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GENETICS_VCF_PHASE2_HPP
#define BOOST_GENETICS_VCF_PHASE2_HPP

#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <cctype>

namespace boost {
namespace genetics {
namespace vcf {

/// @brief Percent encoding/decoding utilities for VCF
namespace percent_encoding {

    /// @brief Encode special characters for VCF INFO/FORMAT fields
    /// Encodes space, semicolon, equals, comma, newline, tab, and percent
    inline std::string encode(const std::string& str) {
        std::ostringstream oss;
        for (std::size_t i = 0; i < str.length(); ++i) {
            char c = str[i];
            if (c == ' ' || c == ';' || c == '=' || c == ',' || 
                c == '\r' || c == '\n' || c == '\t' || c == '%') {
                oss << '%' << "0123456789ABCDEF"[(c >> 4) & 0xF]
                    << "0123456789ABCDEF"[c & 0xF];
            } else {
                oss << c;
            }
        }
        return oss.str();
    }

    /// @brief Decode percent-encoded string
    inline std::string decode(const std::string& str) {
        std::ostringstream oss;
        for (std::size_t i = 0; i < str.length(); ++i) {
            if (str[i] == '%' && i + 2 < str.length()) {
                char hex[3] = {str[i+1], str[i+2], 0};
                int value = 0;
                std::istringstream iss(hex);
                iss >> std::hex >> value;
                oss << static_cast<char>(value);
                i += 2;
            } else {
                oss << str[i];
            }
        }
        return oss.str();
    }

} // namespace percent_encoding

/// @brief Header metadata structures

/// @brief INFO field definition
struct info_meta {
    std::string id;
    std::string number;  // ".", "A", "R", "G", or integer
    std::string type;    // "Integer", "Float", "Flag", "Character", "String"
    std::string description;
    std::map<std::string, std::string> extra; // Additional fields

    info_meta() {}
    
    info_meta(const std::string& id_, const std::string& number_,
              const std::string& type_, const std::string& description_)
        : id(id_), number(number_), type(type_), description(description_) {}
};

/// @brief FORMAT field definition
struct format_meta {
    std::string id;
    std::string number;
    std::string type;
    std::string description;
    std::map<std::string, std::string> extra;

    format_meta() {}
    
    format_meta(const std::string& id_, const std::string& number_,
                const std::string& type_, const std::string& description_)
        : id(id_), number(number_), type(type_), description(description_) {}
};

/// @brief FILTER definition
struct filter_meta {
    std::string id;
    std::string description;

    filter_meta() {}
    
    filter_meta(const std::string& id_, const std::string& description_)
        : id(id_), description(description_) {}
};

/// @brief Contig definition
struct contig_meta {
    std::string id;
    std::size_t length;
    std::map<std::string, std::string> extra;

    contig_meta() : length(0) {}
    
    contig_meta(const std::string& id_, std::size_t length_ = 0)
        : id(id_), length(length_) {}
};

/// @brief ALT allele definition (for structural variants)
struct alt_meta {
    std::string id;
    std::string description;

    alt_meta() {}
    
    alt_meta(const std::string& id_, const std::string& description_)
        : id(id_), description(description_) {}
};

/// @brief Parse structured header line (##KEY=<...>)
/// Returns map of key=value pairs
inline std::map<std::string, std::string> parse_structured_header(const std::string& content) {
    std::map<std::string, std::string> result;
    
    // Remove enclosing < >
    std::string inner = content;
    if (inner.length() >= 2 && inner[0] == '<' && inner[inner.length()-1] == '>') {
        inner = inner.substr(1, inner.length() - 2);
    }
    
    // Parse key=value pairs, respecting quotes
    std::size_t pos = 0;
    while (pos < inner.length()) {
        // Find key
        std::size_t eq_pos = inner.find('=', pos);
        if (eq_pos == std::string::npos) break;
        
        std::string key = inner.substr(pos, eq_pos - pos);
        pos = eq_pos + 1;
        
        // Find value
        std::string value;
        if (pos < inner.length() && inner[pos] == '"') {
            // Quoted value
            ++pos;
            std::size_t end_quote = inner.find('"', pos);
            if (end_quote == std::string::npos) break;
            value = inner.substr(pos, end_quote - pos);
            pos = end_quote + 1;
            
            // Skip comma
            if (pos < inner.length() && inner[pos] == ',') ++pos;
        } else {
            // Unquoted value
            std::size_t comma_pos = inner.find(',', pos);
            if (comma_pos == std::string::npos) {
                value = inner.substr(pos);
                pos = inner.length();
            } else {
                value = inner.substr(pos, comma_pos - pos);
                pos = comma_pos + 1;
            }
        }
        
        result[key] = value;
    }
    
    return result;
}

} // namespace vcf
} // namespace genetics
} // namespace boost

#endif // BOOST_GENETICS_VCF_PHASE2_HPP
