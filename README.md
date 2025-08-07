# Dragonfly-Inspired Microglider 
### Summary
This project was conducted at the Imperial College London Aerial Robotics 
Lab (ARL) under Professor Armanini by Kasper Atkinson (kna25@imperial.ac.uk).
A research group in the ARL was investigating the potential design of a
dragonfly-inspired microglider. The purpose of this project was to look
deeper into the design variables and constraints that would be present in
such a system. A secondary purpose was to determine the aerodynamic
feasibility of controlling the microglider's flight through a C.G.-ramping
mechanism.

### File Structure
Code files are found in the `src` folder, but are organized into key folders
based on which function their script was meant to perform. Here are basic 
overviews of each:
- `design`: broad multi-variable sweeps analyzing their effects on stability and efficiency.
- `general`: multi-purpose functions utilized throughout the whole codebase.
- `linear`: eigenvalue and mode investigation of the flight system.
- `pid`: closed-loop scripts for controlling the model's pitch.
- `point`: basic point-mass implementation of the equations of flight.
- `ramping`: open-loop scripts for controlling the model's pitch.
- `rigid`: preliminary flight model for glider with a tail.
- `stable`: scripts studying the static margin of various configurations.
- `tandem`: close-to-final tandem-wing-based flight model of system.

### Important Functions
Key scripts were used in later versions of the program to eliminate 
repetition and standardize the code across multiple folders.

The `getConstructionVector()` and `getParameterVector()` functions establish
the basic component constructions of the aircraft. Specifically, the 
parameter function contains the aerodynamic coefficients and finds the total
mass and center of gravity of the system. As the codebase progressed, these
were used heavily.

The `xMassForTargetCG()` function enabled the "setting" of a c.g. by
solving for the external added mass position necessary to shift it by a set
amount.

### Warnings
In the final days of the project, the `findTrimCondition()` and
`findEigenValues()` functions were modified for additional outputs. These
changes may not be reflected across all scripts, particularily in earlier
ones. 