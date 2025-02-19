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
- Required Python packages:
  ```
  pysam
  typer
  pandas
  numpy
  ```

### Installing the Package

1. Clone the repository:
```bash
git clone [your-repository-url]
```

2. Install in development mode:
```bash
pip install -e .
```

Or install the required dependencies:
```bash
pip install pysam typer pandas numpy
```



## Usage

The tool is run from the command line with three required arguments:

```bash
python src/strfinder/findstr.py <bam_file> <config_file> <output_file>
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
max_dist = 5
```

### Configuration Files

1. Chromosome Mapping File (CSV):
   - Maps reference names to chromosome names
   - Required columns: `ref`, `chrom`

2. STR Patterns File (CSV):
   - Defines STR patterns to search for
   - Required columns: `chrom`, `name`, `seq`

## Output Format

The tool generates a CSV file with the following columns:
- `STRName`: Name of the STR pattern
- `Motif`: The repeat motif sequence
- `MotifCount`: Number of repeat occurrences

## Examples

1. Basic usage:
```bash
python src/strfinder/findstr.py sample.bam config.toml results.csv
```

2. Example chromosome mapping CSV:
```csv
ref,chrom
chr1,1
chr2,2
```

3. Example STR patterns CSV:
```csv
chrom,name,seq
1,STR1,CAG
2,STR2,AT
```

## Technical Details

The tool:
- Processes BAM files using pysam
- Uses regex for pattern matching
- Combines adjacent motifs within specified distance
- Handles chromosome name mapping
- Supports variable repeat counts

