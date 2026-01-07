//  Copyright (c) 2026 Boost.Genetics
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GENETICS_SEQUENCE_HPP
#define BOOST_GENETICS_SEQUENCE_HPP

#include <string>
#include <algorithm>
#include <stdexcept>

namespace boost {
namespace genetics {

/// @brief Represents a genetic sequence (DNA, RNA, or protein)
class sequence {
public:
    enum class type {
        DNA,
        RNA,
        PROTEIN
    };

    /// @brief Default constructor
    sequence() : data_(), type_(type::DNA) {}

    /// @brief Constructor with sequence data
    /// @param seq The sequence string
    /// @param t The sequence type (default: DNA)
    sequence(const std::string& seq, type t = type::DNA)
        : data_(seq), type_(t) {
        validate();
    }

    /// @brief Get the sequence data
    const std::string& data() const { return data_; }

    /// @brief Get the sequence length
    std::size_t length() const { return data_.length(); }

    /// @brief Get the sequence type
    type get_type() const { return type_; }

    /// @brief Calculate GC content (for DNA/RNA sequences)
    double gc_content() const {
        if (type_ == type::PROTEIN) {
            throw std::runtime_error("GC content not applicable to protein sequences");
        }
        
        std::size_t gc_count = 0;
        for (char c : data_) {
            if (c == 'G' || c == 'C' || c == 'g' || c == 'c') {
                ++gc_count;
            }
        }
        return data_.empty() ? 0.0 : static_cast<double>(gc_count) / data_.length();
    }

    /// @brief Get the reverse complement (for DNA sequences)
    sequence reverse_complement() const {
        if (type_ != type::DNA) {
            throw std::runtime_error("Reverse complement only applicable to DNA sequences");
        }
        
        std::string result = data_;
        std::reverse(result.begin(), result.end());
        
        for (char& c : result) {
            switch (c) {
                case 'A': c = 'T'; break;
                case 'T': c = 'A'; break;
                case 'G': c = 'C'; break;
                case 'C': c = 'G'; break;
                case 'a': c = 't'; break;
                case 't': c = 'a'; break;
                case 'g': c = 'c'; break;
                case 'c': c = 'g'; break;
                default: break; // Leave other characters unchanged
            }
        }
        
        return sequence(result, type::DNA);
    }

private:
    std::string data_;
    type type_;

    void validate() const {
        // Basic validation - can be extended
        if (data_.empty()) {
            return;
        }
        
        const char* valid_chars = nullptr;
        switch (type_) {
            case type::DNA:
                valid_chars = "ATGCNatgcn";
                break;
            case type::RNA:
                valid_chars = "AUGCNaugcn";
                break;
            case type::PROTEIN:
                valid_chars = "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*";
                break;
        }
        
        if (valid_chars) {
            std::string valid_set(valid_chars);
            for (char c : data_) {
                if (valid_set.find(c) == std::string::npos) {
                    throw std::invalid_argument("Invalid character in sequence");
                }
            }
        }
    }
};

} // namespace genetics
} // namespace boost

#endif // BOOST_GENETICS_SEQUENCE_HPP
