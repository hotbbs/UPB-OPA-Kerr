# UPB-OPA-Kerr

This repository contains the numerical simulation codes associated with our manuscript:

**"Photon blockade effect from synergistic optical parametric amplification and driving force in Kerr-medium single-mode cavity"**

Published in: *Acta Physica Sinica*, **74**(16), 164205 (2025). DOI: [10.7498/aps.74.20250712](https://doi.org/10.7498/aps.74.20250712)
English preprint: [arXiv:2601.01819](https://arxiv.org/abs/2601.01819) [quant-ph]

## 📂 Repository Contents

*   **`/figures`**: Contains the final high-resolution image files for all figures presented in the manuscript.
*   **`/codes`**: Contains the MATLAB scripts used for numerical simulations.

## 💻 Computational Environment

The numerical simulations were carried out using MATLAB. The steady state of the Lindblad master equation for the open quantum system is solved using the open-source **Quantum Optics Toolbox for MATLAB** (S. M. Tan, *J. Opt. B*, **1**, 424 (1999)).

*   **Required Software**: MATLAB R2019b or later 
*   **Required Toolbox**: Quantum Optics Toolbox (by Sze M. Tan)

## 🚀 Usage

1. Ensure that the Quantum Optics Toolbox is installed and added to your MATLAB path.
2. (Recommended) Start a parallel pool for the parameter-sweep scripts: `>> parpool('local')`.
3. Run the scripts in the `/codes` directory. The physical parameters (e.g., $F$, $G$, $\Delta_a$, $\phi_a$, $U$, $\kappa$) are pre-defined at the beginning of each script and can be easily modified to reproduce the results.
4. The scripts will output the numerical data and generate the corresponding 2D density plots or 1D line graphs. If `parfor` overhead dominates on your machine (small grids), simply replace `parfor` with `for` — the two-layer loop structure is preserved and results are identical.

## 📝 Citation

If you use these codes or data in your work, please cite our paper:

> **Zhiqiang Zhang**,  
> *Photon blockade effect from synergistic optical parametric amplification
> and driving force in Kerr-medium single-mode cavity*,  
> **Acta Physica Sinica** **74**(16), 164205 (2025).  
> DOI: [10.7498/aps.74.20250712](https://doi.org/10.7498/aps.74.20250712)  
> English version: [arXiv:2601.01819](https://arxiv.org/abs/2601.01819) [quant-ph]

### BibTeX

```bibtex
@article{Zhang2025PhotonBlockade,
  title   = {Photon blockade effect from synergistic optical parametric
             amplification and driving force in Kerr-medium single-mode cavity},
  author  = {Zhang, Zhiqiang},
  journal = {Acta Physica Sinica},
  volume  = {74},
  number  = {16},
  pages   = {164205},
  year    = {2025},
  doi     = {10.7498/aps.74.20250712},
  note    = {English version: arXiv:2601.01819 [quant-ph]}
}
