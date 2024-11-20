# -*- coding: utf-8 -*-
import unittest
from field import Field
import sys
from random import randint
from timeit import timeit


class BenchField:

    def set_up_field(self):
        """Create Bandersnatch base field.

        Test vectors generated using the file `field_test_vectors.sage`.

        """
        try:
            with open("src/tests/field_test_vectors.py", "r") as file:
                exec(file.read(), globals())
        except FileNotFoundError as e:
            raise unittest.SkipTest(
                "The file 'field_test_vectors.py' was not found. Please generate it using `sage field_test_vectors.sage > field_test_vectors.py`.")
        return F, test_vectors  # type: ignore

    def bench_arithmetic(self):
        """Benchmark field multiplication."""
        F, test_vectors = self.set_up_field()
        n_iter = 2000
        naive_mul_time = timeit(
            "tmp=test_vectors['a']*test_vectors['b']", globals=globals(), number=n_iter)
        naive_add_time = timeit(
            "tmp=test_vectors['a']+test_vectors['b']", globals=globals(), number=n_iter)
        naive_div_time = timeit(
            "tmp=test_vectors['a']/test_vectors['b']", globals=globals(), number=n_iter)

        print("Mul: {:.3f}μs".format(naive_mul_time/n_iter * 10**6))
        print("Add: {:.3f}μs".format(naive_add_time/n_iter * 10**6))
        print("Div: {:.3f}μs".format(naive_div_time/n_iter * 10**6))

    def bench_larger_field(self):
        for p in [
            0x949e517f4288a3f1d5402b6230b5aa3c4799c8c4fdf7bad1310b0150620a3d4681a3eed1de7f8b664d7f63dbe27a32944b3620b0a9a7842d8687ae2d4825c3f5,
            0x448a1be8a7514d5097b88d6ad5c92cb12d12753d48a61fc6f0b3594f28353ee18cc723bce5684c67cb3f59fbef65a672679155e23366cc98b24b5043ea2d982ef968e0599aac36394520d2dc7d895cbb0ae0f3f836e695dd20cfcb380e83989738a945c55a91e8845d835bff27ad5c120d6d09210909c350fe0f67c83f8de77b,
            0xbe9078dae23ecdc21dd051f94db7d4ab3b67c013fe087fac73fedf4eb8f8df97e9f857a33f5ef36af8061fc61fd8e32c3ae14bef124e4ac6606f1dccd02d7d9f677fb610c36361dddebc76df0f6093f1e3ae23bab1753fa2f20824155ca1a863cc16f39cb3faf4bded87903398d04cb53f50dff91eaaa627295954f2024b4b26f623e6d1956ac3a8e4518b228d9f72800f129c52efa791ba3304a25c89862d10f129469e0793597d2b6392b183966f8ef1ef997e6eaaeb9c84f8a8a99cc73981fdf18638411dc7d398f80836eda77471c55024b85751a6c8a9a418d32fabbb10a2d1943d13a57e5e8151118344f3d9a7839483c646d52e335c8011b9dd7190fd
        ]:
            print("Benchmark for a prime of size {}".format(p.bit_length()))
            F = Field(p)
            a = F.random()
            b = F.random()
            from time import time
            n_iter = 100000
            t = time()
            for i in range(n_iter):
                _ = a*b
            t_mul = (time()-t)/n_iter

            t = time()
            for i in range(n_iter):
                _ = a+b
            t_add = (time()-t)/n_iter

            t = time()
            for i in range(n_iter):
                _ = a/b
            t_div = (time()-t)/n_iter

            print("\tTime for mul: {:.2f}μs".format(t_mul*10**6))
            print("\tTime for add: {:.2f}μs".format(t_add*10**6))
            print("\tTime for div: {:.2f}μs".format(t_div*10**6))

    def run_all_bench(self):
        for method_name in dir(self):
            if method_name.startswith("bench_"):
                method = getattr(self, method_name)
                if callable(method):
                    print(f"Running: {method_name}")
                    method()
