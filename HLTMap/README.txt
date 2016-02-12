HLTMap.cpp is a program to select sub-areas of primary mirrors in RICH2 in order to select tracks suitable for RICH2 reconstruction.
For each mirror pair, calculates best area to select, using a variable size and shape of binning. Also calculates the reconstruction efficiency and the total amount of events needed to filter through to populate that mirror pair, and relative weighting for each.



To run:
SetupProject Brunel
make
./HLTMap



Basic structure:
-Loop over variables from pre-run tuple to fill a (finely binned) TH3 containing x,y co-ordinates and Phi.

-Loop over a variable merging of the TH3 and extract a reconstruction efficiency for each of the different shaped bins.
-----Reconstruction efficiency calculated by considering the proportion of signal photons in selected bin and the Phi distribution.

-Select highest efficiency value bin from variable merging. Can control how small the binning is (to prevent overbinning), max number of events to filter, max number of bins to merge

-As runs writes to screen a two reference binning values (check against overbinning), the best value within selected criteria, and the best value (regardless of events to filter etc), and saves the relevant plots for each of these.

-At end, writes best values to screen and saves map information in a tuple.


Run time:
- from tuple ~ 2 days
- re-running from pre-run TH3 ~ 2 -> 20 mins (depending on max number of mergings)


NOTE:
- In order to switch between left-hand side mirrors and right hand side mirrors, have to manually comment and uncomment the relevant variables which hold mirror numbers etc.


'Option' variables in code:

 bool fromTuple = false; //If get data from the tuple or from an already run .root file

 int numbMirrPair = 47; //number of mirror pairs to iterate over;default should be 47
 
 int tuplePhotons = treePhoton->GetEntries();

 int rebinXMin = 3; //set the minimum rebinning value that is *saved* - it still runs over the lowest
 
 int rebinYMin = 3; // 
 
 double maxAllowedEvents = 400000000; //max events to filter through that will be considered.
  
 int maxBinMerge = 20; // highest number of bin mergings - default 20

  








