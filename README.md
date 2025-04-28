# inversion_finder

finding inversions in a pangenome graph

## Introduction

I was not able to find a good preëxisting tool for identifying inversions in a pangenome graph: the VCFs made by PGGB and minigraph-cactus do not classify variants but only give the REF and ALT sequences, and although vcfwave purports to do this, it usually segfaults when I run it. However, inversions in a pangenome graph are simply series of nodes that are traversed in opposite directions in different paths, so they ought to be easy to find using just the graph structure without thinking about the actual sequence. This program finds such series of nodes.

A few warnings are in order: I designed this program for a specific task, and it works well for that task in its current form, but it has some big limitations, especially when it comes to finding large inversions. I'm hoping to improve this, but until then, for inversions bigger than 1Mb or so, I recommend using minimap2 in assembly-to-assembly mode and then making a dotplot rather than using this program.

## Installation

This program is written in rust. If you don't have rust installed, you can do this in your user directory with one command:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then, just clone this repository and build the binary:

```bash
git clone https://github.com/esrice/inversion_finder.git
cd inversion_finder
cargo build -r
```

The program will then be in `target/release/inversion_finder`.

## Input preparation

Currently, this program can only be run on a single chromosome at a time, and it requires each chromosome to be in a single scaffold in each assembly. I'm hoping to fix this, but for now, here are the steps to prepare input:

1. Extract the chromosome of interest from each assembly, e.g., with `samtools faidx assembly1.fa chr1`.
2. Make a single fasta file containing all of the sequences to align, where each is named in [PanSN format](https://github.com/pangenome/PanSN-spec), e.g., `assembly1#0#chr1`.
3. Run [PGGB](https://github.com/pangenome/pggb). This may require some futzing with parameters to get the best alignment.

Now, you will have a pangenome in GFA format that you can use as input to this program.

## Running

Now, just run the binary you built:

```bash
inversion_finder pggb_output.gfa name_of_ref_path
```

where `name_of_ref_path` is the name of the assembly you want to use as a reference to compare all the other assemblies to. It can be either the full path name (e.g., `assembly1#0#chr1`) or just the assembly name (e.g., `assembly1`).

The resources required are highly dependent on the structure of the input graph, with graphs with longer inversions present in more assemblies taking longer and requiring more memory. If you run out of memory, try decreasing the `--max-highmem-path-length` parameter. If the program is getting stuck on long inversions (you can see this in the STDERR logging messages), you can try reducing the `--max-path-length` option.

The output is a table of inversions. The first three columns are chromosome, inversion start, and inversion end, in 1-based coordinates of the reference. The rest of the columns are the calls for the non-reference assemblies; a 1 indicates this segment of the assembly is inverted compared to the reference, whereas a 0 indicates it is not.
