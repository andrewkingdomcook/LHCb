Overview:

Program to create tuples with relevant physical quantities (e.g. mass, momentum) of particle objects. 



Basic structure:

- Iterates over particle objects retrieved from LHCb data files. 
- Various selection criteria are applied (e.g. removes particles outside a certain mass range).
- For particles that survive selection criteria, relevant data is written to a tuple for further manipulation.



To run:

Run using the LHCb physics analysis software ‘Davinci’. More information can be found about the Davinci project here:

http://lhcb-release-area.web.cern.ch/LHCb-release-area/DOC/davinci/

Due to the size of the datasets, program is executed using the LHCb distributed computing system, ‘Ganga’.







