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

`rfbzero.py` is a Python package for simulating electrochemical cycling of redox flow batteries. This package includes modules for initial cell setup and electrolyte description, cycling protocol selection, and optional inputs for various capacity degradation mechanisms and active species crossover.
currently allows for zero dimensional modelling of electrochemical cycling techniques and possible capacity fade mechanisms including:
<break this down further?>
- _redox flow cell setup_: easy configuration of flow cell and electrolyte parameters
- _cycling protocol_: quickly define the desired electrochemical cycling protocol
- _capacity fade mechanisms_: include optional degradation and crossover mechanism inherent to the electrolytes and cell


# Introduction
Redox flow batteries (RFBs) are seen as a promising long-duration energy storage technology for grid-scale applications. Zero-dimensional models have previously been developed to understand the electrochemical cycling behaviour of vanadium-based RFBs [@2018_Konig_JPOWERSOURCE; @2018_Murthy_JES; @2018_Pugach_APPENERG; @2019_Lee_JECHEMENERGCONVSTOR], where the dominant capacity fade mechanism involves crossover of active species between negative electrolyte (negolyte) and positive electrolyte (posolyte) reservoirs. Development of second generation electrolyte chemistries, such as redox-active organic molecules (RAOMs) [@2020_Kwabi_CHEMREV], in the past decade requires new models that incorporate properties inherent to novel chemistries. It is often the case that RAOMs are sufficiently bulky so as not to experience appreciable membrane crossover in an RFB, yet unlike vanadium ion electrolytes, they can experience chemical degradation leading to capacity decay. Recent work [@2021_Modak_JES; @2022_Neyhouse_JES] has extended VRFB-based zero-dimensional models to now include the effect of chemical degradation of redox active organics in RFBs.

# Statement of Need
To date, zero-dimensional RFB models have typically been disseminated in the literature via ad hoc non-generalizable code and often written in proprietary software languages. With `rfbzero.py` we hope to provide an open-source Python package that accomplishes electrochemical learning objectives for RFBs, as well as allows for expansion of battery diagnostics via understanding of capacity fade mechanisms observed in the RAOM flow battery community.

The initial cell setup configuration can often dictate the future trend in temporal capacity.

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



