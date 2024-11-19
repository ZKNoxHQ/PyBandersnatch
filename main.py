# -*- coding: utf-8 -*-
import argparse

from curve_bench import BenchCurve
from curve_test import TestCurve
from field_test import TestField


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Run tests or benchmarks for the Curve.")
    parser.add_argument('--test', action='store_true',
                        help="Run tests for the curve")
    parser.add_argument('--bench', action='store_true',
                        help="Run benchmarks for the curve")

    args = parser.parse_args()

    if args.test:
        test_curve = TestCurve()
        test_curve.run_all_test()
        test_field = TestField()
        test_field.run_all_test()
    elif args.bench:
        bench_curve = BenchCurve()
        bench_curve.run_all_bench()
    else:
        print("Please specify either --test or --bench.")
