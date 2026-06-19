# gcutils

`gcutils` is a small Python utility package developed for this project. It provides shared functions for processing Gubbins output, working with homology maps, visualizing gene conversion events on phylogenies, and other genomics tasks used across multiple analysis notebooks.

## Installation

From the repository root:

```bash
pip install -e .
```

## Modules

| Module | Description |
|---|---|
| `gubbinsfuncs.py` | Functions for parsing and processing Gubbins recombination prediction output (GFF, branch statistics CSV, node-labelled phylogeny) |
| `homologymapfuncs.py` | Functions for loading and querying the minimap2 H37Rv self-alignment homology map; maps genomic intervals to paralogous regions |
| `homologymap_graph.py` | Graph-based analysis of paralogous region relationships using NetworkX; constructs paralog networks for PE/PPE and other multigene families |
| `eventtoparalogcomparison.py` | Functions for mapping Gubbins-predicted gene conversion events to their corresponding paralogous donor/acceptor regions |
| `parsimony.py` | Parsimony-based ancestral state reconstruction on phylogenetic trees; used to infer the directionality of gene conversion events |
| `genomeviz_utils.py` | Genome visualization utilities for plotting variants, gene conversion events, and coverage across the H37Rv reference |
| `treeviz_ete_utils.py` | Phylogenetic tree visualization using [ETE3](http://etetoolkit.org/); functions for annotating trees with gene conversion events and lineage metadata |
| `treeviz_mpl.py` | Phylogenetic tree visualization using matplotlib; lightweight alternative to ETE for embedding trees in figure panels |
| `kmer.py` | K-mer based sequence utilities |
| `general.py` | General-purpose helper functions used across modules (I/O, data manipulation) |
