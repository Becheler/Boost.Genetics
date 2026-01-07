//  Copyright (c) 2026 Boost.Genetics
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <genetics/sequence.hpp>
#include <benchmark/benchmark.h>
#include <string>

using namespace boost::genetics;

static void BM_SequenceCreation(benchmark::State& state) {
    const std::string seq_data = "ATCGATCGATCGATCGATCGATCGATCGATCG";
    for (auto _ : state) {
        sequence seq(seq_data, sequence::type::DNA);
        benchmark::DoNotOptimize(seq);
    }
}
BENCHMARK(BM_SequenceCreation);

static void BM_GCContent(benchmark::State& state) {
    sequence seq("ATCGATCGATCGATCGATCGATCGATCGATCG", sequence::type::DNA);
    for (auto _ : state) {
        double gc = seq.gc_content();
        benchmark::DoNotOptimize(gc);
    }
}
BENCHMARK(BM_GCContent);

static void BM_ReverseComplement(benchmark::State& state) {
    sequence seq("ATCGATCGATCGATCGATCGATCGATCGATCG", sequence::type::DNA);
    for (auto _ : state) {
        sequence rc = seq.reverse_complement();
        benchmark::DoNotOptimize(rc);
    }
}
BENCHMARK(BM_ReverseComplement);

static void BM_LongSequenceGCContent(benchmark::State& state) {
    std::string long_seq;
    for (int i = 0; i < state.range(0); ++i) {
        long_seq += "ATCG";
    }
    sequence seq(long_seq, sequence::type::DNA);
    
    for (auto _ : state) {
        double gc = seq.gc_content();
        benchmark::DoNotOptimize(gc);
    }
}
BENCHMARK(BM_LongSequenceGCContent)->Range(8, 8<<10);

BENCHMARK_MAIN();
