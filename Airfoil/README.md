## Project Description
The project simulated the NACA 0012 airfoil (shown in Figure 1) in an outer flow fields of 45 m in the x-direction and 60 m in the y-direction. For use in fluent meshing, the geometry is extruded in the z-direction a distance of 1/10 of the chord length and symmetry conditions are applued to the side boundaries.

The airfoil was simulated under the following operating conditions:
</br> Adhesion angle = 1.55 degrees
</br> Mach number = 0.7
</br> Temperature = 311K
</br> Pressure 101325 Pa
</br> The results are used to compare standard experimental data for verification (shown in Figure 2).The simulation report from Ansys is attached for reference.

![airfoil](https://user-images.githubusercontent.com/110358483/213298355-76254747-07f0-4110-97fe-923202fa9110.png)
Figure 1: Detail of NACA 0012 airfoil  (http://airfoiltools.com/airfoil/details?airfoil=naca0021-il)
</br>

![airfoil_sol](https://user-images.githubusercontent.com/110358483/213297615-df7f1f57-670e-439b-8a9f-337ee8521d01.png)
Figure 2: Final output compare with experimental data
</br>

## About meshing

The meshing in this project is for fluids only. The project was meshed using polyhedra to reduce the number of cells and to improve the quality. A local sizing was add for the 0.0007m gap in Figure 1, and the mesh size of upper and lower layers of the airfoil was set to be 0.015m.

https://user-images.githubusercontent.com/110358483/213306929-d107d56b-7b88-430b-8151-188392efe84a.mp4


