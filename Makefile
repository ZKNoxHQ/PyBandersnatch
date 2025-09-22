.SILENT: clean test benchmark
MAKEFLAGS += --no-print-directory
PY = myenv/bin/python

install:
	python -m venv myenv
	myenv/bin/pip install -r requirements.txt

clean:
	find . -type f -name '*.sage.py' -exec rm -f {} +
	find . -type d -name '__pycache__' -exec rm -rf {} +

gen_test_vec:
	# test vectors for Bandersnatch 
	sage sage/bandersnatch_field.sage > tests/vectors/bandersnatch_field.py
	sage sage/bandersnatch_montgomery.sage > tests/vectors/bandersnatch_montgomery.py
	sage sage/bandersnatch_edwards.sage > tests/vectors/bandersnatch_edwards.py
	# test vectors for Ed25519
	sage sage/ed25519_field.sage > tests/vectors/ed25519_field.py
	make clean

test:
	@if [ -z "$(TEST)" ]; then \
		$(PY) -m unittest discover -s tests; \
	else \
		$(PY) -m unittest tests.$(TEST); \
	fi
	make clean

benchmark:
	@if [ -z "$(BENCH)" ]; then \
		$(PY) -m unittest discover -s bench -q; \
	else \
		$(PY) -m unittest bench.$(BENCH) -q; \
	fi
	make clean
