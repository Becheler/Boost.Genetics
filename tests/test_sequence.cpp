//  Copyright (c) 2026 Boost.Genetics
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <genetics/sequence.hpp>

using namespace boost::genetics;

TEST_CASE("sequence basic functionality", "[sequence]") {
    SECTION("default constructor") {
        sequence seq;
        REQUIRE(seq.data().empty());
        REQUIRE(seq.length() == 0);
    }

    SECTION("constructor with data") {
        sequence seq("ATCG", sequence::type::DNA);
        REQUIRE(seq.data() == "ATCG");
        REQUIRE(seq.length() == 4);
        REQUIRE(seq.get_type() == sequence::type::DNA);
    }

    SECTION("invalid DNA sequence throws") {
        REQUIRE_THROWS_AS(sequence("ATCGX", sequence::type::DNA), std::invalid_argument);
    }
}

TEST_CASE("sequence GC content", "[sequence]") {
    SECTION("all GC") {
        sequence seq("GCGC", sequence::type::DNA);
        REQUIRE(seq.gc_content() == Approx(1.0));
    }

    SECTION("no GC") {
        sequence seq("ATAT", sequence::type::DNA);
        REQUIRE(seq.gc_content() == Approx(0.0));
    }

    SECTION("50% GC") {
        sequence seq("ATGC", sequence::type::DNA);
        REQUIRE(seq.gc_content() == Approx(0.5));
    }

    SECTION("empty sequence") {
        sequence seq("", sequence::type::DNA);
        REQUIRE(seq.gc_content() == Approx(0.0));
    }

    SECTION("GC content with lowercase") {
        sequence seq("atgc", sequence::type::DNA);
        REQUIRE(seq.gc_content() == Approx(0.5));
    }

    SECTION("protein sequence throws") {
        sequence seq("ACDEFGH", sequence::type::PROTEIN);
        REQUIRE_THROWS_AS(seq.gc_content(), std::runtime_error);
    }
}

TEST_CASE("sequence reverse complement", "[sequence]") {
    SECTION("simple sequence") {
        sequence seq("ATCG", sequence::type::DNA);
        sequence rc = seq.reverse_complement();
        REQUIRE(rc.data() == "CGAT");
    }

    SECTION("palindromic sequence") {
        sequence seq("GAATTC", sequence::type::DNA);
        sequence rc = seq.reverse_complement();
        REQUIRE(rc.data() == "GAATTC");
    }

    SECTION("lowercase sequence") {
        sequence seq("atcg", sequence::type::DNA);
        sequence rc = seq.reverse_complement();
        REQUIRE(rc.data() == "cgat");
    }

    SECTION("RNA sequence throws") {
        sequence seq("AUCG", sequence::type::RNA);
        REQUIRE_THROWS_AS(seq.reverse_complement(), std::runtime_error);
    }

    SECTION("protein sequence throws") {
        sequence seq("ACDEFGH", sequence::type::PROTEIN);
        REQUIRE_THROWS_AS(seq.reverse_complement(), std::runtime_error);
    }
}
