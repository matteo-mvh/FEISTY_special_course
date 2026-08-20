# FEISTY Fishing & Carbon Model

> **Project type:** DTU University Research Project / Special Course  
> **Field:** Marine ecosystem modelling · biological carbon pump · fisheries · numerical modelling  
> **Tools:** FEISTY · MATLAB / numerical modelling

## Overview

This project investigates how **fishing alters marine ecosystem structure and the biological carbon pump** using the FEISTY size-structured marine ecosystem model.

The work combines local ecosystem experiments with global simulations to examine how changes in fish biomass, size structure and functional groups can influence the transport and long-term storage of carbon in the ocean.

The repository was developed as part of a **DTU special course in marine ecological modelling**.

## Research Questions

The project focuses on questions such as:

- How does fishing change fish biomass and food-web structure?
- Do different fishing strategies produce different effects on carbon export?
- How do these effects vary between shelf seas, slopes and the open ocean?
- What is the potential global impact of fishing on biologically mediated carbon transport and storage?

## Project Structure

The repository is divided into two complementary components:

```text
FEISTY_special_course/
├── Local/
├── Global/
└── README.md
```

### Local simulations

The `Local/` component investigates idealised ecosystems representing different marine environments:

- **Shelf sea (~75 m)** – productive coastal ecosystem with strong demersal and forage-fish components
- **Continental slope (~1500 m)** – intermediate-depth ecosystem
- **Open ocean (~3000 m)** – deep, lower-productivity environment

These simulations are used to understand mechanisms before moving to global-scale analysis.

### Global simulations

The `Global/` component applies FEISTY at global scale and combines ecosystem simulations with environmental forcing and ocean transport information.

The objective is to evaluate whether fishing-induced changes in marine biomass can alter the amount and distribution of carbon transported into the ocean interior.

## Fishing Scenarios

The model is used to compare an unfished baseline with alternative fishing strategies, including fishing pressure targeting different ecological groups such as:

- demersal fish
- forage fish
- large pelagic fish
- size-selective fishing scenarios

Changes are evaluated relative to the unfished ecosystem.

## Main Outputs

The simulations provide quantities including:

- equilibrium fish biomass
- biomass distribution among functional groups and size classes
- food-web responses to fishing
- fishing-induced carbon injection
- differences relative to an unfished baseline
- spatial patterns in global biomass and carbon fluxes
- estimates of carbon sequestration associated with different scenarios

## Modelling Approach

The project follows a two-scale workflow:

1. Explore ecological mechanisms using simplified local ecosystems.
2. Compare ecosystem responses across depth and productivity regimes.
3. Introduce alternative fishing scenarios.
4. Quantify resulting changes in biomass and carbon transport.
5. Extend the analysis to global FEISTY simulations.
6. Combine biological carbon injection with ocean transport to investigate longer-term carbon storage.

## Why This Project Matters

Fish influence marine carbon cycling through feeding, respiration, mortality, vertical movement and the production of sinking organic material.

Fishing therefore has the potential to affect not only fish stocks but also ecosystem-mediated carbon transport.

This project explores that connection using a mechanistic ecosystem model and provides experience in:

- marine ecosystem modelling
- size-structured ecological models
- fisheries impact analysis
- biological carbon cycling
- scenario analysis
- large-scale environmental modelling
- scientific data analysis and visualisation

## Status

Completed university special-course project.

The repository contains research code developed for scientific exploration rather than a production software package.
