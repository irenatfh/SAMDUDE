# SAMdude
Universal denoiser for genomic sequencing

Irena Fischer-Hwang ihwang@stanford.edu

These scripts denoise aligned and sorted SAM files by changing nucleotide bases and updating corresponding quality scores.

## Installation

SAMdude requires no additional software beyond Python 3+ and NumPy.

## Usage

The input SAM files must be sorted. Denoising affects only the read and quality score strings; all headers, identifiers and other fields (flags, CIGAR strings, etc.) will be kept intact.

To denoise a sam file, use the following command:

```
python3 run_denoiser.py noisy_file.sam denoised_file.sam
```
