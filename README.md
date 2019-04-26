# PyrInt-model

## Description
Network model of Hodgkin-Huxley type pyramidal cells and interneurons as used in Ter Wal & Tiesinga, Front. Comput. Neurosci., 2017 (see https://doi.org/10.3389/fncom.2017.00006). The model generates 1 or 2 circuits containing  pyramidal cells and/or interneurons at specified connections probabilities and strengths. The model comes with many input options: static currents, pulses, oscillatory inputs, correlated and uncorrelated noise, or any combination of those, and inputs can be controlled per region. The code is made to loop through different input conditions and/or different connection parameters (see DetermineLoops.m for all options).

## Instructions
To run the model:
1. adjust settings.m. Choose the desired model settings and the path for saving data and figures in the settings.m file;
2. run RunScript.m. You can change what to analyse and what to plot in this file. 

Note: the current model is set up to produce figure 2 of the 2017 paper, i.e. with 2 regions and looping through different oscillation frequencies for region 1 and connection strengths between the regions.

## Disclaimer
The model itself should generalize fairly well to other questions than the ones studied in the paper, and the model comes with more options for input settings than used in the paper. However, keep in mind that the analysis and plotting functions were written for the specific goals of the paper, and will not always work with other configurations of the model. 

The code is still in serious need of commenting and documenting, feel free to help me with this by pointing out the parts that are particularly foggy. 

If you have any comments or suggestions, feel free to raise an issue or contact me at m.j.terwal@bham.ac.uk.
