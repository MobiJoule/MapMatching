# MapMatching

This repository contains code for matching trajectories (of geographic coordinates) to the transport network. It is based on code originally developed by
Jan-Henrik Haunert and Axel Forsch at https://gitlab.igg.uni-bonn.de/graphlibrary/gl-mapmatching, which is based on the methodology described in Haunert and Budig (2012)
(https://doi.org/10.1145/2424321.2424402) and Newson and Krumm (2009) (https://doi.org/10.1145/1653771.1653818).

## Installation
1. Clone the gl-full repository from https://gitlab.igg.uni-bonn.de/graphlibrary/gl-full.git
2. Replace the remote repository for the gl-mapmatching module with this one
3. Run using the main method in matching/main/MatchingMain.java. Use the -h flag to see instructions and arguments.


## About
This repository extends Haunert and Forsch's work to:
* Match to both non-directed and directed graphs.
* Use the latest version of GeoTools (currently 33.1) to support reading/writing shapefiles and geopackages.
* Integrate partitioning (to manage memory) and multi-threaded computation (to decrease runtime),
  allowing large trajectory datsets to be processed at once.
* Output all results to a single geopackage. Outputs are written to the disk incrementally (after each partition) to
  preserve results in case a computation is terminated prematurely.
* Eliminate zero-length segments from output,
* Improve matching accuracy at the start point by using a larger matching weight
  (i.e., penalty for being further from a link) at the start point relative to
  other points along the trajectory.
* Allow modifying link weights based on their attributes. This enables the matching algorithm to favour or avoid
  certain routes when multiple options align with the trajectory (e.g., cycling paths parallel to roads).
  The adjustments are specified in an optional text file read in as an input using the flag "-adj".
  In this text file, each line specifies using the format "attribute = criteria #.#". For example, the following
  modifications might be used for matching cycling trajectories to avoid nearby pedestrian paths:

```aiignore
highway = pedestrian 1.5 
highway = footway 1.5
```

# Work in progress:
* Create flag to specify directed / non-directed graph
* Incorporate timestamps in output (for later estimating speed)
* Make independant and eliminate reliance on old "gl-*" libraries
* Update or remove Ma**r**chingMain.java