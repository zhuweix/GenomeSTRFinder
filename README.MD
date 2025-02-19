# STR Analysis Tool

A command-line tool for analyzing Short Tandem Repeats (STRs) from BAM files. This tool processes aligned sequences and identifies STR patterns based on user-defined configurations.

## Features

- Analyze BAM files for STR patterns
- Configure custom motifs and parameters via TOML
- Support for chromosome name mapping
- Flexible pattern matching with customizable distance parameters
- Output results in CSV format

## Installation

### Prerequisites

- Python 3.11 or later
- Required packages will be automatically installed through pip

### Installing for Development

1. Clone the repository:
```bash
git clone https://github.com/zhuweix/GenomeSTRFinder.git
```

2. Install development dependencies and the package in editable mode:
```bash
cd GenomeSTRFinder
pip install .
```

This will automatically install all required dependencies:
- Core dependencies: pysam, typer, pandas, numpy



## Usage

After installation, you can run the tool using the command-line interface:

```bash
strfinder <bam_file> <config_file> <output_file>
```

Arguments:
- `bam_file`: Path to input BAM file
- `config_file`: Path to TOML configuration file
- `output_file`: Path for output CSV file

## Configuration

The tool uses a TOML configuration file for settings. Example configuration:

```toml
[Data]
chrom_file = "path/to/chromosome_mapping.csv"
str_file = "path/to/str_patterns.csv"

[Parameter]
min_count = 3
max_dist = 150
```

### Configuration Files

1. Chromosome Mapping File (CSV):
   - Maps reference names to chromosome names
   - Required columns: `ref`, `chrom`

2. STR Patterns File (CSV):
   - Defines STR patterns to search for
   - Required columns: `chrom`, `name`, `seq`

   The config/STR_NIST.csv file contains a list of STR motif from the [NIST STRBase database](https://strbase.nist.gov/).

### Parameters

min_count: Minimum number of repeat occurrences to consider a pattern
max_dist: Maximum distance between adjacent motifs to combine

## Output Format

The tool generates a CSV file with the following columns:
- `STRName`: Name of the STR pattern
- `Motif`: The repeat motif sequence
- `MotifCount`: Number of repeat occurrences (N1;N2...;Nk if mutliple motifs are within max_dist distance, set in the configuration)

## Examples

1. Basic usage:
```bash
strfinder sample.bam config.toml results.csv 
```

- results.csv.gz is also supported for compressed output


2. Example STR patterns CSV:
```csv
STRName,Motif,MotifCount
D1S1677,TTCC,4;4
```

## Technical Details

The tool:
- Processes BAM files using pysam
- Uses regex for pattern matching
- Combines adjacent motifs within specified distance
- Handles chromosome name mapping
- Supports variable repeat counts


## License

MIT License
