# Loblaw Bio - miraclib cell-count analysis.
#
#   make setup      create .venv and install the pinned dependencies
#   make pipeline   Part 1: build cell_counts.db  ->  Parts 2-4: write outputs/
#   make dashboard  serve the interactive dashboard on http://localhost:8501
#   make test       run the test suite

PYTHON ?= python3
VENV   := .venv
PY     := $(VENV)/bin/python
PORT   ?= 8501
# GitHub Codespaces serves the dashboard through its port-forwarding proxy, whose requests
# fail Streamlit's CORS / XSRF origin checks, so the local server turns them off. A hosted
# deployment (Streamlit Community Cloud) does not use these flags and keeps the protections.
STREAMLIT_FLAGS := --server.headless true --browser.gatherUsageStats false \
                   --server.enableCORS false --server.enableXsrfProtection false

.DEFAULT_GOAL := help
.PHONY: help setup pipeline dashboard test clean

help: ## List the available targets
	@grep -E '^[a-z]+:.*## ' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*## "}; {printf "  make %-10s %s\n", $$1, $$2}'

setup: ## Create .venv and install the pinned dependencies
	$(PYTHON) -m venv $(VENV)
	$(PY) -m pip install --quiet --upgrade pip
	$(PY) -m pip install --quiet -r requirements.txt
	@echo "Dependencies installed in $(VENV)/ - next: make pipeline"

# Any target that needs Python runs setup first if .venv does not exist yet.
$(PY):
	@$(MAKE) --no-print-directory setup

pipeline: | $(PY) ## Build cell_counts.db from cell-count.csv, then write the Part 2-4 outputs to outputs/
	$(PY) load_data.py
	$(PY) run_analysis.py

dashboard: | $(PY) ## Serve the dashboard (Parts 2-4) at http://localhost:8501 (change with PORT=...)
	@test -f cell_counts.db || $(PY) load_data.py
	$(PY) -m streamlit run dashboard.py --server.port $(PORT) $(STREAMLIT_FLAGS)

test: | $(PY) ## Run the test suite
	$(PY) -m pytest

clean: ## Remove the database and Python caches (outputs/ and .venv are kept)
	rm -f cell_counts.db cell_counts.db.tmp
	rm -rf .pytest_cache
	find . -path ./$(VENV) -prune -o -name __pycache__ -type d -prune -exec rm -rf {} +
