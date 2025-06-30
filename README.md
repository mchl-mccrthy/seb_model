# A simplified energy-balance model for ice melt below debris

These scripts and functions implement a simplified energy-balance model for simulating ice melt below debris at the surfaces of debris-covered glaciers. The model solves the one-dimensional heat equation to conduct heat through the debris layer, following Reid and Brock (2010), where the boundary condition at the debris surface is a simplified debris-surface energy balance (cf. Oerlemans, 2001, for debris-free ice melt), which contains two free parameters to be calibrated. The boundary condition at the ice surface is the temperature of melting ice. The model is described mathematically [here](/docs/seb_model_description.pdf). The model is written in MATLAB and was tested in MATLAB 2022a.

Author: Michael McCarthy (michael.mccarthy@wsl.ch)

## Running the model
The model is set up and run with the script [src/run_seb.m](src/run_seb.m). 

## Input data 
Some synthetic input data are stored in the folder [/inputs](/inputs) as an example. 

## Output data
Running the model will produce output data in the MATLAB terminal, which can be saved according to the user's preferences.

## References
Oerlemans, J. (2001). Glaciers and climate change. CRC Press.

Reid, T. D., & Brock, B. W. (2010). An energy-balance model for debris-covered glaciers including heat conduction through the debris layer. Journal of Glaciology, 56(199), 903-916.

## Referencing the model
The model can be referenced e.g., in APA format as follows.

McCarthy, M. (2025). A simplified energy-balance model for ice melt below debris (Version 1.0.0). GitHub. [https://github.com/mchl-mccrthy/seb_model](https://github.com/mchl-mccrthy/seb_model)
