.PHONY: init test demo

init:
	uv venv
	uv pip install .

dev:
	uv venv
	uv pip install .[dev]
	@echo "Virtual environment created/updated. Activate it by running: source .venv/bin/activate"

test:
	uv run python -m unittest test.test_basic

demo:
	uv run python test/test_demo.py
