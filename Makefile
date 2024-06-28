simulate:
	sh ./sim.sh

docs:
	quarto render downs.ipynb --to html --output-dir ./docs

.PHONY: docs
