# CQ_tools

CQ_tools is a comprehensive bioinformatics pipeline designed for calculating and visualizing **CQ (Chromosome Quotient)** values from sequencing data. It supports genome alignment, automated window-based coverage calculation, and result visualization through both static images and interactive web reports.

## Features

- **Automated Alignment**: Integrated BWA MEM alignment and Samtools processing.
- **Windowed Analysis**: Flexible genomic window generation and coverage calculation using `bedtools`.
- **CQ Calculation**: Efficient CPM (Counts Per Million) normalization and CQ ratio calculation.
- **Unified Logging**: Centralized logging system for tracking progress and debugging.
- **Rich Visualization**:
  - High-resolution static plots (PNG/PDF).
  - Interactive HTML reports with hover details and zooming capabilities.

## Installation

### Prerequisites

- Python >= 3.11
- [BWA](http://bio-bwa.sourceforge.net/)
- [Samtools](http://www.htslib.org/)
- [bedtools](https://bedtools.readthedocs.io/)

### Setup

We recommend using `uv` for dependency management:

```bash
git clone https://github.com/yourusername/CQ_tools.git
cd CQ_tools
uv sync
```

Or install via `pip`:

```bash
pip install pandas matplotlib plotly pysam typer
```

## Usage

CQ_tools provides three main commands: `align`, `cq`, and `plot`.

### 1. Alignment (`align`)

Align reads to a reference genome and generate sorted BAM files.

```bash
python main.py align 
    --fasta reference.fasta 
    --pair-1 sample_R1.fastq.gz 
    --pair-2 sample_R2.fastq.gz 
    --output-dir ./results 
    --cpu 8
```

### 2. CQ Analysis (`cq`)

Calculate coverage, normalize data, and identify specific genomic regions.

```bash
python main.py cq 
    --fasta reference.fasta 
    --f-bam female_sorted.bam 
    --m-bam male_sorted.bam 
    --output-dir ./cq_results 
    --cq-value 0.3 
    --threshold 10 
    --parallel 2
```

### 3. Visualization (`plot`)

Generate distribution plots from existing CQ results.

```bash
python main.py plot 
    --cq-result cq_results/F_M_CQ.filter.tsv 
    --chrom-length cq_results/chromosome_length.tsv 
    --output distribution.png 
    --html
```

#### Example Output Plot:
![CQ Distribution Plot](CQPlot/CQ_test_plot.png)

## Output Files

- `CQ.log`: Detailed execution log.
- `chromosome_windows.tsv`: Generated genomic windows (BED format).
- `F_M_CQ.filter.tsv`: Final filtered CQ results.
- `CQ_distribution.png`: Static genomic distribution plot.
- `CQ_distribution.html`: Interactive web-based report.

## License

MIT License
