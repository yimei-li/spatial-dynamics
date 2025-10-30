# DIP–virus coinfection simulations

**Key insight:**  
Coinfection requires spatial overlap between DIPs and virions.  
**If DIPs jump too far from virions, no coinfection occurs, and because DIPs cannot replicate, the infection collapses entirely.**  
The model reproduces this stochastic threshold — most runs fail, but some yield the experimental pattern.


Validation of our stochastic model against experimental data from  
**Baltes et al. (2017)** — *Inhibition of infection spread by co-transmitted defective interfering particles* (*PLOS ONE*, 12:e0184029).


---


### Figure 2
\[
\text{Model validation in the absence of IFN.}
\]


- Solid lines: experimental data.  
- Dashed lines: mean of stochastic simulations.  
- Grey shading: 95% range of simulation outcomes.


**Key insight:**  
Coinfection requires spatial overlap between DIPs and virions.  
**If DIPs jump too far from virions, no coinfection occurs, and because DIPs cannot replicate, the infection collapses entirely.**  
The model reproduces this stochastic threshold — most runs fail, but some yield the experimental pattern.


---


### Run
```bash
python run_sim.py --config configs/baltes_validation.yaml