# -*- coding: utf-8 -*-
import argparse
from bench.curve_bench import BenchCurve
from bench.curve_edwards_bench import BenchCurveEdwards
from bench.field_bench import BenchField
from curve import Curve
from field import Field
from key_exchange import KeyExchange
from tests.curve_edwards_test import TestCurveEdwards
from tests.curve_test import TestCurve
from tests.field_test import TestField
from tests.key_exchange_test import TestKeyExchange
import sys
import os

# Add the project root to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def run_specific_function(args, classes):
    for cls in classes:
        if hasattr(cls, args.function):
            func = getattr(cls, args.function)
            func()
            return
    print(f"Function {args.function} not found in any provided classes.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Run tests or benchmarks for Bxndxrsnxtch.")
    parser.add_argument('--test', action='store_true',
                        help="Run tests for the curve")
    parser.add_argument('--bench', action='store_true',
                        help="Run benchmarks for the curve")
    parser.add_argument('--function', type=str,
                        help='Specify a function to test/bench')

    args = parser.parse_args()

# Classes to consider
    test_classes = [TestField(), TestCurve(
    ), TestCurveEdwards(), TestKeyExchange()]
    bench_classes = [BenchField(), BenchCurve(), BenchCurveEdwards()]

    if args.function:
        run_specific_function(args, test_classes + bench_classes)
    else:
        # Run all tests or benchmarks if no function is specified
        if args.test:
            for cls in test_classes:
                cls.run_all_test()
        elif args.bench:
            for cls in bench_classes:
                cls.run_all_bench()
