//  Copyright (c) 2026 Boost.Genetics
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <genetics/genetics.hpp>
#include <iostream>

int main() {
    using namespace boost::genetics;

    // Create a SNP variant
    variant snp("chr1", 12345, "A", "G");
    
    std::cout << "Variant at " << snp.chromosome() << ":" << snp.position() << std::endl;
    std::cout << "REF: " << snp.ref() << " -> ALT: " << snp.alt() << std::endl;
    std::cout << "Is transition: " << (snp.is_transition() ? "Yes" : "No") << std::endl;
    std::cout << "Is transversion: " << (snp.is_transversion() ? "Yes" : "No") << std::endl;
    
    // Create an insertion
    variant ins("chr2", 67890, "A", "ATG");
    std::cout << "\nInsertion at " << ins.chromosome() << ":" << ins.position() << std::endl;
    std::cout << "REF: " << ins.ref() << " -> ALT: " << ins.alt() << std::endl;
    
    return 0;
}
