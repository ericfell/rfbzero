---
title: 'RFBzero: A Python package for zero-dimensional simulation of redox flow battery cycling'
tags:
  - Python
  - electrochemistry
  - redox flow batteries
  - capacity fade
authors:
  - name: Eric M. Fell
    email: efell@g.harvard.edu
    affiliation: 1
    orcid: 0000-0003-2046-1480
  - name: Jeremy A. Fell
    affiliation: 2
    orcid: 0009-0001-0968-8564
  - name: Michael J. Aziz
    email: maziz@harvard.edu
    affiliation: 1
    orcid: 0000-0001-9657-9456
affiliations: 
  - name: Harvard John A. Paulson School of Engineering and Applied Sciences, 29 Oxford Street, Cambridge, MA, 02138, USA
    index: 1
  - name: Simon Fraser University, 8888 University Dr W, Burnaby, BC V5A 1S6, Canada
    index: 2
  
date: 28 August 2023
bibliography: paper.bib
---

# Overview

`rfbzero` is a Python package for simulating electrochemical cycling of redox flow batteries. This package includes modules for initial cell setup, cycling protocol selection, and optional inputs for various chemical degradation mechanisms and active species crossover.
currently allows for zero dimensional modelling of electrochemical cycling techniques and possible capacity fade mechanisms including:
<break this down further?>
- _CC_: constant current cell cycling
- _CCCV_: constant current followed by constant voltage cell cycling
- _degradation_: first or second order capacity degradation mechanisms
- _crossover_: crossover of redox-active species

# Introduction
placeholder

# Statement of Need
The initial cell setup configuration can often dictate the future trend in temporal capacity.
Previous work [@2021_Modak_JES] has extended VRFB-based zero-dimensional models to include the effect of chemical degradation of redox active organics in RFBs.

# Current `rfbzero.py` Functionality


## Cell design
The initial cell configuration allows for standard variables which a researcher would normally adapt. Cell specific parameters such as active area, resistance, reservoir volumes of the capacity limiting side (CLS) and non-capacity limiting side (NCLS) are considered in `redox_flow_cell.py`. The initial concentration of redox active species (oxidized and/or reduced), and choice of symmetric cell (identical species both sides, 0 V OCV) or full cell (different species on each side, OCV > 0 V) are also considered.

## Cycling protocol

## Degradation mechanisms

## Crossover mechanism

# An Example of the `rfbzero.py` API
The documention includes...

```python
# an example
something = 6.7
print(something)
```


# Acknowledgements
We thank Prof. David Kwabi, Thomas George, and Jordan Sosa for constructive feedback and testing.

# References



