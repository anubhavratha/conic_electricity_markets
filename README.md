# Moving from linear to conic markets for electricity: Code and Data

## About:
This repository contains the code and simulation data for the paper, [Moving from linear to conic markets for electricity](https://doi.org/10.1016/j.ejor.2022.12.025), published in European Journal of Operational Research (EJOR). A paywall-free final manuscript is also available on [ArXiv](https://arxiv.org/abs/2103.12122).

If you use this code or parts of it, please cite the paper as follows:
> @article{Ratha2023, title = {Moving from linear to conic markets for electricity},
author = {Anubhav Ratha and Pierre Pinson and Hélène {Le Cadre} and Ana Virag and Jalal Kazempour},
journal = {European Journal of Operational Research}, volume = {309}, number = {2}, pages = {762-783}, year = {2023}, issn = {0377-2217}, doi = {https://doi.org/10.1016/j.ejor.2022.12.025}}


## Setup:

The optimization models were implemented in [Julia](https://julialang.org) (v.1.6) using [JuMP](https://github.com/JuliaOpt/JuMP.jl) modeling language for mathematical optimization. The conic electricity market clearing problem is solved using the [Mosek](https://www.mosek.com) solver, which requires a license (free of charge for academic use). Refer to JuMP documentation 

The non-convex models are solved using the [Ipopt](https://ipoptjl.readthedocs.io/en/latest/ipopt.html) solver, while the [Mosek](https://www.mosek.com) solver is used for the convex second-order cone programming models. The Mosek solver needs to be licensed (free of charge for academic use).  Please refer to [JuMP documentation](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) on how to set up solvers. The packages used and their versions are provided in the `Project.toml` file. To activate the packages in ```Project.toml```, open a terminal, clone the project using ```git clone```, ```cd``` to the project directory and run the following code block:
```
$ julia 
julia> ]
pkg> activate .
pkg> instantiate
```

```
julia> include("main.jl")
```

A chance-constrained conic electricity market clearing problem solved as a central planner problem. Further details on the model parameters and the reformulation of the chance constraints can be found in the manuscript.The figure below shows the energy prices at various nodes of the 24-bus electricity system considered in the case study.

<img width="381" alt="Prices" src="https://user-images.githubusercontent.com/19344128/172192699-5064902d-8b1b-4002-87a1-5203a66cf906.png">


