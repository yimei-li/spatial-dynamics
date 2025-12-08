# Representative Simulation Results for Fig 2

We ran the stochastic spatial simulation many times (40 runs) and got a wide range of outcomes. Since it's stochastic, each run can look pretty different depending on where virions and DIPs happen to land and spread.

We picked 6 representative snapshots (`represent1.png` ~ `represent6.png`) to show this variability.

## IMP: Why the results vary so much

The big thing is: **DIPs can't replicate on their own**!
They need to coinfect a cell together with a functional virus to hijack its replication machinery. So what happens early on matters a lot:
- If initial DIPs happen to land near viruses and infect the same cells, you get coinfection spreading and yellow patches show up
- If DIPs land somewhere far from viruses and never meet them, those DIPs just sit there doing nothing - they can't spread, can't coinfect, basically they're dead weight
- Sometimes DIPs get lucky and catch up with the virus wave later, sometimes they don't


This is why some simulation runs show lots of yellow (coinfected) patches 
while others are mostly red (virus-only) with 
barely any DIP (green) involvement.


## Files
- `simulation_25_hours_*.png` - All individual simulations, I extracted them at t=25 hours
- `represent1.png` - `represent6.png` - Selected representative outcomes for our figure 2  
- `all.png` - Combined view
- `main.go` - The simulation code
