//  Copyright (c) 2026 Boost.Genetics
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <genetics/genetics.hpp>
#include <iostream>

int main() {
    using namespace boost::genetics;

    // Create a DNA sequence
    sequence dna("ATCGATCG", sequence::type::DNA);
    
    std::cout << "Sequence: " << dna.data() << std::endl;
    std::cout << "Length: " << dna.length() << std::endl;
    std::cout << "GC Content: " << dna.gc_content() * 100 << "%" << std::endl;
    
    // Get reverse complement
    sequence rc = dna.reverse_complement();
    std::cout << "Reverse Complement: " << rc.data() << std::endl;
    
    return 0;
}
