PYTHON ?= python3

.PHONY: test demo doctor
test:
	$(PYTHON) -m unittest discover -s tests -v

demo:
	$(PYTHON) -m protein_stability demo --out runs/demo

doctor:
	$(PYTHON) -m protein_stability doctor
