ARG PYTHON_IMAGE=python:3.12-slim
FROM ${PYTHON_IMAGE}
WORKDIR /app
COPY protein_stability /app/protein_stability
COPY examples /app/examples
COPY tests /app/tests
# No pip download and no proprietary software in the image.
CMD ["python", "-m", "unittest", "discover", "-s", "tests", "-v"]
