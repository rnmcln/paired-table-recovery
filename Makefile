PY = python3

.PHONY: all screen analyse validate tables figures test clean

all: validate analyse tables figures

screen:
	$(PY) scripts/01_screen.py

analyse:
	$(PY) scripts/02_analyse.py

validate:
	$(PY) scripts/03_validate.py

figures:
	$(PY) scripts/04_figures.py

tables:
	$(PY) scripts/05_tables.py

test:
	$(PY) -m pytest -q tests

clean:
	rm -rf results/*.csv results/*.txt figures
