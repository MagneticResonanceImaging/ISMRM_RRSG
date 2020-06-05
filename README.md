# ISMRM 2019 Reproducible research study group challenge

### Description

Implementation of [ISMRM Reproducibility Challenges](https://ismrm.github.io/rrsg/) using the MRI reconstruction framework [MRIReco.jl](https://travis-ci.org/github/MagneticResonanceImaging/MRIReco.jl). 

### Challenge 1

Goal of [challenge 1](https://blog.ismrm.org/2019/04/02/ismrm-reproducible-research-study-group-2019-reproduce-a-seminal-paper-initiative/) is to reproduce the results from the paper [1]. The julia scripts and results can be found in the subfolder *challenge1*. To run the code one first has to download [julia](https://julialang.org/downloads/). Then, follow the [instructions](https://magneticresonanceimaging.github.io/MRIReco.jl/latest/#Installation-1) to install MRIReco.jl. Finally run one of the scripts by entering *include("reco_....jl")* into the Julia terminal. To do so, one either needs to startup julia in the folder *challenge1* or use the [shell mode](https://docs.julialang.org/en/v1/stdlib/REPL/#man-shell-mode-1) of julia to *cd* into that folder. 

#### Results
The following figures show the reconstruction results generated with *MRIReco.jl*.

![Figure 4](challenge1/Fig4.png?raw=true "Figure 4")
![Figure 5](challenge1/Fig5.png?raw=true "Figure 5")

### References

[1] Pruessmann, K. P.; Weiger, M.; Boernert, P. and Boesiger, P. Advances in sensitivity encoding with arbitrary k-space trajectories. Magn Reson Med 46: 638-651 (2001)

### Contact

Tobias Knopp [tobias . knopp (at) tuhh . de] and Mirco Grosser [mirco . grosser (at) tuhh . de]

[Institute for Biomedical Imaging](tuhh.de/ibi),   
University Medical Center Hamburg Eppendorf and Hamburg University of Technology, Germany.


