# pairdown
A CLI tool to downsample (pare down, get it?) **paired-end** FASTQs written in Rust.

I couldn't find one that preserved read pairing in two FASTQ files for some reason, though I'm sure one exists. If you do find one, don't tell me, otherwise I'll feel like I wasted my time doing this.


## Installation

Make sure you have Rust installed, then run:

```
cargo install --git https://github.com/quoctran98/pairdown.git
```

## Usage

Run:

```
pairdown read1.fq read2.fq 12000000 -o downsampled_ -s 1
```

Output:

```
Downsampled 12000000 from 14887030 reads to downsampled_read1.fq and downsampled_read2.fq in 110 seconds.
```


This command takes in two paired-end FASTQ files (`read1.fq` and  `read2.fq`), randomly samples (without replacement, using a random seed `1`) 12,000,000 reads, then writes them to two output FASTQ files (`downsampled_read1.fq` and  `downsampled_read2.fq`).

## Options

```
-o, --prefix      Output file prefix (default: downsampled_<reads_kept>_)
-s, --seed        Random seed for reproducibility (default: 17)
```
