## Project Description
In this project, we simulated the heat transfer of the heat sink in a container. We also compared the case with and without radiation and the results are shown in the figure below. The simulation report from Ansys is attached for reference.

<img src="https://user-images.githubusercontent.com/110358483/214164063-71c9f21b-3d37-407b-a700-cca2db983ab6.png" width=40% height=40%> <img src="https://user-images.githubusercontent.com/110358483/214164089-c782a173-bcff-44ed-b754-136be3df4923.png" width=50% height=50%>

Figure 1: heat transfer of heat sink with (left) and without(right) radiation

<img src="https://user-images.githubusercontent.com/110358483/214166498-1be23c27-92d4-484e-b8a6-89b7440a32b0.png" width=40% height=40%> <img src="https://user-images.githubusercontent.com/110358483/214166575-479dc581-cc7f-42ed-8b04-18ee4f003f1b.png" width=40% height=40%>

Figure 2: heat transfer of container with (left) and without(right) radiation

The Figure 1 shows that the temperature of the heat sink is 14 degrees lower than the temperature with the radiation module turned off. And the radiation of the system has a significant effect on the the simulation output of the heat sink container. Therefore, opening the radiation module in this simulation is necessary.

## About meshing 
Local sizing were added to the entire heat sink, which significantly increased the number of cells. In this project, a polyhexagonal mesh was used to meet the cell limits, the details of which are given below. Compared to hexcore, poly-hexcore has fewer elements and, in general, uses the least amount of RAM.



https://user-images.githubusercontent.com/110358483/214229988-3921b8a8-b1a9-4ab0-aed3-4e453f247d2c.mp4

