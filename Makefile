.SILENT: clean test benchmark
MAKEFLAGS += --no-print-directory

install:
	pip install -r requirements.txt

clean:
	find . -type f -name '*.sage.py' -exec rm -f {} +
	find . -type d -name '__pycache__' -exec rm -rf {} +

gen_test_vec:
	sage sage/field.sage > tests/vectors/field.py
	sage sage/montgomery.sage > tests/vectors/montgomery.py
	sage sage/edwards.sage > tests/vectors/edwards.py
	make clean

test:
	@if [ -z "$(TEST)" ]; then \
		python -m unittest discover -s tests; \
	else \
		python -m unittest tests.$(TEST); \
	fi
	make clean

benchmark:
	@if [ -z "$(BENCH)" ]; then \
		python -m unittest discover -s bench -q; \
	else \
		python -m unittest bench.$(BENCH) -q; \
	fi
	make clean
