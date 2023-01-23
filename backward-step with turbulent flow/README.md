## Project Description

This project simulates the motion pattern of turbulent flow through pipes of different diameters. The figure below shows the cross-sectional view of the pipe and the fluid direction for this project. The project was considered to have a constant velocity input of 25 m/s. The actual experimental velocity input and shear stress output were recorded for further analysis.

<img src="https://user-images.githubusercontent.com/110358483/214159427-e9778a5a-30ea-4c21-a57a-b22fc7a997bf.png" width=70% height=70%>

The results are considered accurate due to the mass flow rate at the inlet and the mass flow rate at the outlet shown below, and the residuals of the solution meet the desired value. The simulation report from Ansys is attached for reference.


                  Mass Flow Rate               [kg/s]
                ----------------- --------------------
                           inlet          0.015897868
                          outlet         -0.015897867
                ---------------- --------------------
                             Net        8.1239903e-10
                             
                             
This project also compares the shear stresses along the y-direction with the constants (purple line) and the experimental inputs (green line) to the experimental output. 

<img src="https://user-images.githubusercontent.com/110358483/214156392-4fe90aba-f7cd-4ade-a5e6-3717071bdffb.png" width=70% height=70%>

## About meshing
The transition between two different diameters is considered to be a critical part of this simulation. Therefore, we added local sizing for meshing in both Y-direction and X-direction for that section as below.

