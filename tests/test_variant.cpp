//  Copyright (c) 2026 Boost.Genetics
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <genetics/variant.hpp>

using namespace boost::genetics;

TEST_CASE("variant basic functionality", "[variant]") {
    SECTION("default constructor") {
        variant v;
        REQUIRE(v.chromosome().empty());
        REQUIRE(v.position() == 0);
        REQUIRE(v.ref().empty());
        REQUIRE(v.alt().empty());
    }

    SECTION("constructor with data") {
        variant v("chr1", 12345, "A", "G");
        REQUIRE(v.chromosome() == "chr1");
        REQUIRE(v.position() == 12345);
        REQUIRE(v.ref() == "A");
        REQUIRE(v.alt() == "G");
    }
}

TEST_CASE("variant type detection", "[variant]") {
    SECTION("SNP") {
        variant v("chr1", 100, "A", "G");
        REQUIRE(v.get_type() == variant::type::SNP);
    }

    SECTION("insertion") {
        variant v("chr1", 100, "A", "ATG");
        REQUIRE(v.get_type() == variant::type::INSERTION);
    }

    SECTION("deletion") {
        variant v("chr1", 100, "ATG", "A");
        REQUIRE(v.get_type() == variant::type::DELETION);
    }

    SECTION("indel - same length but not SNP") {
        variant v("chr1", 100, "AT", "GC");
        REQUIRE(v.get_type() == variant::type::INDEL);
    }
}

TEST_CASE("variant transition and transversion", "[variant]") {
    SECTION("A to G is transition") {
        variant v("chr1", 100, "A", "G");
        REQUIRE(v.is_transition());
        REQUIRE_FALSE(v.is_transversion());
    }

    SECTION("G to A is transition") {
        variant v("chr1", 100, "G", "A");
        REQUIRE(v.is_transition());
        REQUIRE_FALSE(v.is_transversion());
    }

    SECTION("C to T is transition") {
        variant v("chr1", 100, "C", "T");
        REQUIRE(v.is_transition());
        REQUIRE_FALSE(v.is_transversion());
    }

    SECTION("T to C is transition") {
        variant v("chr1", 100, "T", "C");
        REQUIRE(v.is_transition());
        REQUIRE_FALSE(v.is_transversion());
    }

    SECTION("A to C is transversion") {
        variant v("chr1", 100, "A", "C");
        REQUIRE_FALSE(v.is_transition());
        REQUIRE(v.is_transversion());
    }

    SECTION("A to T is transversion") {
        variant v("chr1", 100, "A", "T");
        REQUIRE_FALSE(v.is_transition());
        REQUIRE(v.is_transversion());
    }

    SECTION("insertion is neither") {
        variant v("chr1", 100, "A", "ATG");
        REQUIRE_FALSE(v.is_transition());
        REQUIRE_FALSE(v.is_transversion());
    }

    SECTION("deletion is neither") {
        variant v("chr1", 100, "ATG", "A");
        REQUIRE_FALSE(v.is_transition());
        REQUIRE_FALSE(v.is_transversion());
    }
}
