# PPG07 Analysis: $\pi^0/\eta$ Transverse Single-Spin Asymmetries in $p^\uparrow + p$ collisions

This is the analysis code of the PPG07 analysis ([Conf Note](https://www.sphenix.bnl.gov/PublicResults/sPH-CONF-COLDQCD-2025-02)).

The code is divided into three sections:

  - [Run Selection](runSelection): Scripts for extracting the good run list from the SQL database.
  - [Trigger Bit Selection](triggerBitSelection): Pre-select events for which the photon and MBD scaled trigger bits are fired
  - [Main Analysis](mainAnalysis):
    - Extract relevant cluster information from the pre-selected triggered events and store them.
    - Extract candidate photon pair/meson candidates + QA plots, after passing all cuts.