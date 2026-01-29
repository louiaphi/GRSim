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
- Be sick for a week 😔//sry for using Emojis, im not ChatGPT frfr
- 28.1.2026 Made function to construct Camera Tetrad and tested it... didn't work 😔 again
- 19.1.2026 Made Heatmap Visualization of H field Strength. Inside the ring singularity r ≈ 0 -> circle singularity (Note 1)

## ToDo

- check for constraints in geodesic integration [x]
- chech for constraints in all dimensions [ ]
- create local tedrad [ ]
- account for colissions [ ]
- render image [ ]
- account for redshift [ ]

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
