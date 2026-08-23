# Supported Cluster Generator

A Python package for generating supported mono- and multimetallic cluster structures.

The current version focuses on topological growth from an existing supported cluster. It reads structures with ASE, identifies cluster metals, builds a NetworkX graph of the metal framework, finds growth motifs and generates new one-centre structures in either surface-parallel or 3D directions.

Final version aim: A topology-aware, ML-driven structure-discovery framework for supported multimetal clusters, designed to efficiently explore the enormous structural and compositional spaces encountered in heterogeneous catalysis.

## Current functionality

- Read one file, multiple files, directories or glob patterns
- Identify user-defined cluster metals
- Calculate metal-metal distances using MIC
- Build a NetworkX metal framework
- Analyse neighbours, coordination, edges and triangles
- Identify one-, two- and three-centre growth motifs
- Generate surface-parallel or 3D growth directions
- Generate one- two- and three centre growth
- Filter uphysical candidates
- Write new ASE structures
- Run from a command-line interface

## Installation

From the repository root:

```bash
python -m pip install -e ".[dev]"
```

## Example

Surface-parallel growth:

```bash
cluster-gen grow inputs/trimer.xyz --metals Co Mn Ru --add Ru --geometry surface
```

3D growth:

```bash
cluster-gen grow inputs/trimer.xyz --metals Co Mn Ru --add Ru --geometry 3d --directions 16
```
