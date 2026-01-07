//  Copyright (c) 2026 Boost.Genetics
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <genetics/variant.hpp>
#include <benchmark/benchmark.h>

using namespace boost::genetics;

static void BM_VariantCreation(benchmark::State& state) {
    for (auto _ : state) {
        variant v("chr1", 12345, "A", "G");
        benchmark::DoNotOptimize(v);
    }
}
BENCHMARK(BM_VariantCreation);

static void BM_VariantIsTransition(benchmark::State& state) {
    variant v("chr1", 12345, "A", "G");
    for (auto _ : state) {
        bool is_trans = v.is_transition();
        benchmark::DoNotOptimize(is_trans);
    }
}
BENCHMARK(BM_VariantIsTransition);

static void BM_VariantIsTransversion(benchmark::State& state) {
    variant v("chr1", 12345, "A", "C");
    for (auto _ : state) {
        bool is_transv = v.is_transversion();
        benchmark::DoNotOptimize(is_transv);
    }
}
BENCHMARK(BM_VariantIsTransversion);

BENCHMARK_MAIN();
