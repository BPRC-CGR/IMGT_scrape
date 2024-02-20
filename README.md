# IMGT Sequence Scraper README

## Overview

This script automates scraping of immunogenetics (IG) and T-cell receptor (TCR) VDJ segment sequences from the International ImMunoGeneTics Information System (IMGT) for specified species. Which is specified in the **--help**. It supports fetching sequences in various open reading frame (ORF) analysis types and consolidates them into a single sequence library. Additionally, it offers a cleanup **--cleanup** function to remove the individual fasta files by deleting downloaded these after processing.

## Features

- Species-specific scraping: Fetches IG and TCR segment sequences for a specified species.
- ORF frame selection: Offers options for different ORF frame analysis types (all, in-frame, in-frame-gaps).
- Library creation: Allows for the creation of a consolidated sequence library from the fetched files.
- Cleanup: Option to automatically delete fetched files after processing.

## Requirements

- Python 3.11 or higher (Configured in **env.yaml** file)
- Required Python libraries: requests, bs4 (BeautifulSoup), Bio (Biopython) (Configured in **env.yaml** file)
- Internet connection for IMGT website access

## Environment Setup

This script requires specific Python packages. An environment.yaml (**env.yaml**) file is provided to easily create a conda environment with all necessary dependencies.

## Creating the Conda Environment

You must have conda installed! If not installed: Visit Conda Installation Guide for instructions.

Create the environment, this is straightforward: Navigate to the project directory, probably called IMGT_scrape with the terminal and run the following command to create the conda environment from the env.yml file:

```bash
conda env create -f envs/env.yaml
```

Activate the environment: Once the environment is created, activate it with:

``` bash
conda activate IMGT
```

## Usage

The script provides a command-line interface with several options. Here is a description of the available arguments:

``` txt
Required Arguments
-S, --species: Specify the species to scrape for (e.g., "Homo sapiens"). Script handles capitalization automatically.
-T, --type: Type of sequence to scrape: TCR (T-cell receptor) or IG (Immunoglobulin).

Optional Arguments

-O, --output: Output directory for the results. Defaults to a directory based on the species name.
-f, --frame-selection: ORF frame analysis type (all, in-frame, in-frame-gaps).
--create-library: Consolidate fetched sequences into a single library file.
--cleanup: Delete fetched files after processing.
```

## Example Command

```python3 IMGT_scrape.py -S "Homo sapiens" -T TCR --create-library --cleanup```

This command fetches TCR sequences for Homo sapiens, creates a library from them, and cleans up the workspace.

## Logging

The script uses custom logging with colored output for clear and structured execution feedback. It supports DEBUG, INFO, WARNING, ERROR, and CRITICAL levels for detailed progress tracking and troubleshooting.

## Contributions

Contributions are welcome. Please fork the repository, make your changes, and submit a pull request.