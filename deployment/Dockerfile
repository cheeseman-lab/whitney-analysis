FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends git && rm -rf /var/lib/apt/lists/*

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy visualization code
COPY brieflow/visualization/ ./visualization/
COPY brieflow/workflow/lib/shared/file_utils.py ./workflow/lib/shared/file_utils.py

# Create __init__.py files for the workflow package path
RUN touch ./workflow/__init__.py ./workflow/lib/__init__.py ./workflow/lib/shared/__init__.py

# Initialize a git repo so the Analysis Overview page doesn't crash
RUN git config --global user.email "deploy@whitney-viz" && git config --global user.name "deploy" && git init && git add -A && git commit -m "deployment" --allow-empty

# Copy baked-in data (~2.4 GB: cluster, aggregate TSVs, eval files)
COPY data/ ./data/

# Environment variables
# BRIEFLOW_OUTPUT_PATH points to a merged view set up by entrypoint
ENV BRIEFLOW_OUTPUT_PATH=/app/merged_data/brieflow_output/
ENV CONFIG_PATH=/app/data/config/config.yml
ENV SCREEN_PATH=/app/data/screen.yaml

# Streamlit config
RUN mkdir -p /root/.streamlit
RUN echo '[server]\nheadless = true\nport = 8080\naddress = "0.0.0.0"\nenableCORS = false\nenableXsrfProtection = false\n' > /root/.streamlit/config.toml

COPY entrypoint.sh .
RUN chmod +x entrypoint.sh

EXPOSE 8080

ENTRYPOINT ["./entrypoint.sh"]
