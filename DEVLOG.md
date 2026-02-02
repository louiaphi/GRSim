# Plan

- Raytracing um ein Schwarzes Loch
- Calc metric Tensor
- Calc Cristophels
- Integrate Geodesic Equation numerically
- Check for hit
- Übertragung auf GPU
- Render rückwärtslaufen lassen

# Change of Plans

- No GPU Acceleration
- Simulate Light one after another
- Simulate Redshift

# Devlog

- 13.1.2026 Calc Metric Tensor
- 15.1.2026 Dependencies to Calculate Metric Tensor
- 16.1.2026 Function to calculate derivatives of Metric Tensor
- 18.1.2026 Runge Kutta RK4 Integration finished
- 20.1.2026 Checked for Constraints and started debugging
- Be sick for a week 😔
<!--- sry for using Emojis, im not ChatGPT frfr--->
- 28.1.2026 Made function to construct Camera Tetrad and tested it... didn't work 😔 again
- 29.1.2026 Made Heatmap Visualization of H field Strength. Inside the ring singularity r ≈ 0 -> circle singularity (Note 1)
- 1.2 Started Camera Debugging (If you include "making it worse" in the definition of debugging)

## ToDo

- check for constraints in geodesic integration [x]
- chech for constraints in all dimensions [ ]
- create local tedrad [x]
- account for colissions [ ]
- render image [ ]
- account for redshift [ ]
- make steps inverse proportional to H so that stronger gravity = smaller step size

## Debugg Steps:

- metric tensor [x]
- christoffel symbols [x]
- del metric tensor [x]
- RK4 [x]
- Visualize Light Stepping [ ]

## Notes

### Note 1: Heatmap

M = 8, a = 4, z = 0, y = 0, x is the variable. 
At x = 0: r = 0, H = NaN 
At 0<x<4: r = 0/very small/NaN, H = very big/ NaN 
At x = 4: r = 0; H = NaN At x > 4: r increases like expected, H rapidly droppes of (like expected)

<br />
<br />
<br />
<!--- brrrrrrrrrr --->

## Interesting / Funny Pictures / Problems encountered while debugging (that will take up way to much space because i such at md)
### Struggles with the GPU:
#### Longhole  <br />
<img width="880" height="493" alt="Screenshot 2026-01-29 110921" src="https://github.com/user-attachments/assets/4f3de198-9565-43c7-92ad-42d482a59599" />  <br />
#### Diagonalhole  <br />
<img width="354" height="372" alt="Screenshot 2026-01-29 104931" src="https://github.com/user-attachments/assets/ff1046a8-a820-466f-8ebc-7a258f0335ab" />  <br />
#### Many Worlds? Na, Many Holes!  <br />
<img width="735" height="739" alt="Screenshot 2026-01-29 161322" src="https://github.com/user-attachments/assets/51b0aff4-8950-4f40-aa4d-0f9dbfae8afb" />  <br />



