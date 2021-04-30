Things to add: 

1. Flux sync up between levels 
done 2. Simplified embedded boundary (base memos version) 
done 3. 3D solution 
done 4. Material properties and phase change 
5. 

SW height and momentum dummy solvers, height continuity for given melt velocity implemented 
Energy redistribution 
Failcheck 2D solver (does not work presently) 
Possible non-uniform grid refinement (only refine in y-direction for instance). Reading AMReX forum seems possible in Cpp impl. but when 
I try in fortran interface (implement grid refinement ratio as vector, unique ratio for each level and spatial direction), error message comes up and run is aborted. 





Notes on strategies: 

In current implementation of surface boundary, the flux is hindered ONLY on surface faces. This means that heat still diffuses in 'vacuum' region above surface. 
Idea to implement an integer multifab with 3 components (for each spatial direction), nodal for the spatial dir. Value 1 on material faces and value 0 on vacuum faces. Multiply face flux with integer multifab value. Thus there is no diffusion in vacuum region. 
This would make it possible to use highest refinement only on molten parts of the surface, and low res across the rest of the interface. 
Possible problem: When interpolating data to fine grid from coarse upon addition of fine grid, energy is 'left hanging' in the vacuum region just above the surface. could reset temperature when defining new level from coarse. 

use algorithm to 

 




