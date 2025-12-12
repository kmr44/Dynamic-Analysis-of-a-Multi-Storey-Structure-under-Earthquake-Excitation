# Dynamic Analysis of a Multi-Storey Structure under Earthquake Excitation

This project performs the dynamic analysis of a 5-DOF shear-frame structure subjected to the 1940 El-Centro earthquake record. The system is modelled using lumped masses and inter-storey stiffness values, and solved using modal analysis to study vibration behaviour, damping influence, and response amplification across floors.

---

## 🔹 Key Features
- **5-DOF lumped-mass shear building model**
- **Modal analysis** to compute:
  - Natural frequencies
  - Mode shapes
  - Modal participation factors
- **Damping incorporated** via modal damping ratios
- **MATLAB implementation** of:
  - Mass and stiffness matrix formulation
  - Eigenvalue-based uncoupling
  - Duhamel’s integral for time-domain response
- **Ground motion input:** 1940 El-Centro earthquake dataset
- **Outputs include:**
  - Displacement, velocity, and acceleration response of each storey
  - Inter-storey drift and vibration amplification trends
  - Comparison with published seismic benchmark behaviour

---

## 📁 Repository Structure

- **Code/**
  - `main.m` — MATLAB script

- **Data/**
  - `elcentro.csv` — Ground acceleration dataset

- **CAD/**
  - `5_story_structure_model.SLDPRT` — CAD model of the structure

- **Report/**
  - `Dynamic_Analysis_of_Multi_Story_Structure_Report.pdf` — Project report

- **README.md** — Project documentation


