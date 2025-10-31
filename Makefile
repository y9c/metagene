.PHONY: help install test demo build clean lint format docs

help:  ## Show this help message
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'


install:  ## Install the package and dependencies
	uv sync

install-dev:  ## Install with development dependencies
	uv sync --extra dev

test:  ## Run tests
	uv run python -m pytest test/ -v

test-cov:  ## Run tests with coverage
	uv run python -m pytest test/ --cov=markdup --cov-report=html --cov-report=term

demo: ## Run demo script
	@uv run python test/test_demo.py

build:  ## Build the package
	uv build

clean:  ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf htmlcov/
	rm -rf .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

lint:  ## Run linting
	uv run ruff check .
	uv run ruff format --check .

format:  ## Format code
	uv run ruff format .
	uv run ruff check --fix .

docs:  ## Build documentation (placeholder)
	@echo "Documentation is in the docs/ directory"
	@echo "To build with mkdocs: mkdocs build"
	@echo "To serve with mkdocs: mkdocs serve"

check: lint test  ## Run all checks

all: clean install-dev test build  ## Run full pipeline
