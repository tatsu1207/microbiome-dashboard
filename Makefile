# ============================================================================
# MicrobiomeDash — Makefile
# ============================================================================

.PHONY: help setup run db-reset clean test lint

ENV_NAME = microbiome
PYTHON = conda run -n $(ENV_NAME) python
UVICORN = conda run -n $(ENV_NAME) uvicorn

help: ## Show this help message
	@echo "MicrobiomeDash — Available commands:"
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'
	@echo ""

setup: ## Run the setup script
	chmod +x setup.sh
	./setup.sh

run: ## Start the application
	chmod +x run.sh
	./run.sh

db-reset: ## Reset the database (WARNING: deletes all data records, not files)
	@echo "⚠️  This will delete the database and recreate it empty."
	@read -p "Are you sure? (y/N): " confirm; \
	if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
		rm -f microbiome.db; \
		$(PYTHON) -c "from app.db.database import init_db; init_db()"; \
		echo "✓ Database reset."; \
	else \
		echo "Cancelled."; \
	fi

clean: ## Remove temporary files, exports, and __pycache__
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	rm -rf data/exports/*
	@echo "✓ Cleaned temporary files."

clean-all: ## Remove ALL data (uploads, datasets, exports, database)
	@echo "⚠️  This will delete ALL data including uploads and datasets."
	@read -p "Are you sure? (y/N): " confirm; \
	if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
		rm -rf data/uploads/* data/datasets/* data/combined/* data/exports/*; \
		rm -f microbiome.db; \
		echo "✓ All data removed."; \
	else \
		echo "Cancelled."; \
	fi

test: ## Run tests
	$(PYTHON) -m pytest tests/ -v

lint: ## Run linting
	$(PYTHON) -m flake8 app/ --max-line-length 120 --ignore E501,W503

check-env: ## Verify environment and dependencies
	@echo "Checking environment..."
	@conda run -n $(ENV_NAME) python -c "import fastapi, dash, plotly, sqlalchemy, pandas, numpy; print('✓ Python packages OK')"
	@conda run -n $(ENV_NAME) Rscript -e 'library(dada2); message("✓ R packages OK")'
	@echo "✓ Environment check passed."
