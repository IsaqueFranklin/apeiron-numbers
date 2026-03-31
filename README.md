# Apeiron: Atomic-Scale Entropy via High-Entropy Alloys (HEAs)

**Apeiron** is a research framework designed to investigate the feasibility of using the intrinsic chemical disorder of **High-Entropy Alloys (HEAs)** and **Ultra-High-Entropy Alloys (UHEAs)** as a physical source of entropy for cryptographic hardware. By mapping the stochastic fluctuations of the **Potential Energy Landscape (PES)** at the atomic level, this project develops a methodology for generating high-quality **True Random Numbers (TRNG)**.

---

## 🔬 Scientific Overview

Traditional hardware TRNGs often rely on thermal noise or ring oscillators. **Apeiron** shifts this paradigm to the material's atomic architecture. In HEAs (e.g., **CoCrFeNi**), the "chemical chaos" created by multiple principal elements results in a rugged energy landscape. We treat this lattice as a natural stochastic medium where carrier transport and atomic transitions become inherently unpredictable due to local neighborhood variations.

### Key Objectives
* **Atomic Modeling:** Constructing FCC/BCC supercells for HEAs and UHEAs using the **Atomic Simulation Environment (ASE)**.
* **Physical Grounding:** Integrating real-world thermodynamic data (enthalpies of mixing, electronegativities) sourced from the **AFLOW** database to ensure simulations reflect physical reality.
* **Stochastic Dynamics:** Utilizing **Kinetic Monte Carlo (kMC)** to simulate transition processes governed by local energy barriers.
* **Cryptographic Validation:** Subjecting the extracted bitstreams to the **NIST SP 800-22** statistical test suite to verify randomness quality.

---

## 🛠 Methodology

The project is structured into three primary computational layers:

### 1. Construction of the Stochastic Medium
We generate crystal lattices where elements are distributed to mimic a solid solution state. By increasing the number of principal elements (moving from HEAs to UHEAs), we analyze how configurational entropy scales with the quality of the generated randomness.

### 2. Data-Driven Energy Mapping
Instead of using idealized mathematical models, Apeiron utilizes **Data Mining via the AFLOW API**. For every atomic site, the local environment’s energy is calculated based on:
* Interatomic interaction parameters.
* Diffusion barriers specific to the alloy's composition.
* Lattice distortion effects.

### 3. Entropy Extraction & Post-Processing
The core entropy source is derived from the **waiting times** and **transition probabilities** of carriers jumping between sites. Because every site has a unique chemical neighborhood, these transitions occur at erratic intervals.
* **Hopping Mechanism:** Simulated via kMC, where transitions are thermally activated and barrier-dependent.
* **Digitization:** The erratic timing of these physical events is converted into binary sequences.
* **Bias Correction:** We implement a **Von Neumann Extractor** to eliminate any statistical bias (unbalance between 0s and 1s) inherent in the physical process.

---

## 📊 Statistical Validation

To prove that "Atomic Chaos" is viable for cybersecurity, all generated data is benchmarked against the **National Institute of Standards and Technology (NIST) SP 800-22** suite. This includes tests for:
* Frequency and Block Frequency.
* Cumulative Sums.
* Longest Run of Ones.
* Non-overlapping Template Matching.

---

## 🚀 Tech Stack

* **ASE (Atomic Simulation Environment):** For lattice construction and manipulation.
* **AFLOW API:** For retrieving thermodynamic and crystallographic properties.
* **Python:** Main engine for kMC implementation and data orchestration.
* **NIST SP 800-22:** Statistical validation framework.

---

> **Etymology:** *Apeiron* (Greek: ἄπειρον) refers to the "infinite" or "indefinite" primordial substance from which all things arise. In this project, it represents the boundless complexity of the atomic landscape as a source of pure randomness.
