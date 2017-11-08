# thin_disk

* [hipparcos_data_cut.py](https://github.com/jtbuch/thin_disk/blob/master/hipparcos_data_cut.py) is my implementation of the 
kinematic cuts on Hipparcos data (without extinction correction) as described in Eric's paper.

* [PJSolver_dark_grid_solve.py](https://github.com/jtbuch/thin_disk/blob/master/PJSolver_dark_grid_solve.py) solves the 
Poisson-Jeans (PJ) equation, estimates the error through sampling, and calculates the chi-square value over a non-rectangular 
grid for plotting contours.
   
  [PJSolver_dark_grid_solve.ipynb](http://nbviewer.jupyter.org/github/jtbuch/thin_disk/blob/master/PJSolver_dark_grid_solve.ipynb) is the same code but in a Jupyter script.
  
* [PJSolver_no_disk_v2.ipynb](http://nbviewer.jupyter.org/github/jtbuch/thin_disk/blob/master/PJSolver_no_disk_v2.ipynb) implements the PJ equation for a standard Bahcall model in the absence of a thin disk.
