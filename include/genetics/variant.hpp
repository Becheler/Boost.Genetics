//  Copyright (c) 2026 Boost.Genetics
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GENETICS_VARIANT_HPP
#define BOOST_GENETICS_VARIANT_HPP

#include <string>
#include <cstddef>

namespace boost {
namespace genetics {

/// @brief Represents a genetic variant (SNP, insertion, deletion, etc.)
class variant {
public:
    enum class type {
        SNP,           // Single Nucleotide Polymorphism
        INSERTION,
        DELETION,
        INDEL
    };

    /// @brief Default constructor
    variant() : chromosome_(), position_(0), ref_(), alt_(), type_(type::SNP) {}

    /// @brief Constructor with variant details
    /// @param chr Chromosome identifier
    /// @param pos Position on the chromosome
    /// @param ref Reference allele
    /// @param alt Alternative allele
    variant(const std::string& chr, std::size_t pos, 
            const std::string& ref, const std::string& alt)
        : chromosome_(chr), position_(pos), ref_(ref), alt_(alt) {
        determine_type();
    }

    /// @brief Get the chromosome
    const std::string& chromosome() const { return chromosome_; }

    /// @brief Get the position
    std::size_t position() const { return position_; }

    /// @brief Get the reference allele
    const std::string& ref() const { return ref_; }

    /// @brief Get the alternative allele
    const std::string& alt() const { return alt_; }

    /// @brief Get the variant type
    type get_type() const { return type_; }

    /// @brief Check if this is a transition (A<->G or C<->T)
    bool is_transition() const {
        if (type_ != type::SNP || ref_.length() != 1 || alt_.length() != 1) {
            return false;
        }
        
        char r = ref_[0];
        char a = alt_[0];
        
        return (r == 'A' && a == 'G') || (r == 'G' && a == 'A') ||
               (r == 'C' && a == 'T') || (r == 'T' && a == 'C');
    }

    /// @brief Check if this is a transversion
    bool is_transversion() const {
        return type_ == type::SNP && ref_.length() == 1 && 
               alt_.length() == 1 && !is_transition();
    }

private:
    std::string chromosome_;
    std::size_t position_;
    std::string ref_;
    std::string alt_;
    type type_;

    void determine_type() {
        if (ref_.length() == alt_.length() && ref_.length() == 1) {
            type_ = type::SNP;
        } else if (ref_.length() < alt_.length()) {
            type_ = type::INSERTION;
        } else if (ref_.length() > alt_.length()) {
            type_ = type::DELETION;
        } else {
            type_ = type::INDEL;
        }
    }
};

} // namespace genetics
} // namespace boost

#endif // BOOST_GENETICS_VARIANT_HPP
