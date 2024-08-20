# mbds-in-bds
Accompanying Python code to the paper [Single-cell mutational burden distributions in birth-death processes](https://arxiv.org/abs/2309.06355), coauthored with Dudley Stark and Weini Huang. Each file is briefly described below, with details in the paper.

## dd.py (with folder 'wijs')

Simulates a birth-death process and plots the single-cell division distribution (DD) averaged over many realisations along with the predicted theoretical distribution. It unpickles stored predicted values from the folder 'wijs' as needed.

## mbd.py

Simulates a birth-death process and plots the single-cell mutational burden distribution (MBD) averaged over many realisations along with the predicted distribution generated with the division distribution (DD).

## pop.py

Plots the expected population size conditioned on survival for a discrete time birth-death process and the linear approximation found in the paper.

## mbd-diff-mut-dists.py

Plots the single-cell mutational burden distribution (MBD) averaged over many realisations, where, depending on which lines are commented out, different mutational distributions can be applied (e.g. geometric or uniform, as done in the paper) instead of Poisson.
