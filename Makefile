.PHONY: init test

init:
	peotry install --no-dev

dev:
	peotry install
	peotry shell

test:
	python -m unittest test.test_basic
