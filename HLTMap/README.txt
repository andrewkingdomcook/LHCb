Overview:
The process of the alignment of certain of the LHCb detector components is slow, because only a small proportion of the data used is suitable for the alignment process. This program calculates a method to pre-select data, so that the time for running the alignment process is greatly reduced.

Specifically, the LHCb ‘RICH’ detector (a basic overview of RICH detectors can be found here https://en.wikipedia.org/wiki/Ring-imaging_Cherenkov_detector) has a considerable number of ‘pairs’ of mirrors, which need to be aligned with each other, in order for the detector to work properly. 

The selection process is based solely on the position of the associated photons when they strike the first set of mirrors in the RICH detector. The program HLTMap.cpp calculates what areas should be used in order to optimise the pre-selection of data.



