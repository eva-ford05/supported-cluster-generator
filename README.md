# Supported Cluster Generator

A Python package for generating supported mono- and multimetallic cluster structures.

The current version focuses on topological growth from an existing supported cluster. It reads structures with ASE, identifies cluster metals, builds a NetworkX graph of the metal framework, finds growth motifs and generates new one-centre structures in either surface-parallel or 3D directions.

## Current functionality

- Read one file, multiple files, directories or glob patterns
- Identify user-defined cluster metals
- Calculate metal-metal distances using MIC
- Build a NetworkX metal framework
- Analyse neighbours, coordination, edges and triangles
- Identify one-, two- and three-centre growth motifs
- Generate surface-parallel or 3D growth directions
- Generate one-centre growth candidates
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

## Planned development

1. Finish one-centre growth
2. Add two-centre edge growth
3. Add three-centre face growth
4. Add clash and support filters
5. Add structural deduplication
6. Add recursive growth to a target size
7. Add bare-surface monomer site search
8. Add a general ASE calculator interface
9. Add MACE relaxation and ranking
10. Add Winterbottom-supported cluster generation
11. Add multimetal composition enumeration
