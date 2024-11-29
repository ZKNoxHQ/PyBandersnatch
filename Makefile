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
	python -m unittest discover -s tests
	make clean

benchmark:
	python -m unittest discover -s bench
	make clean
