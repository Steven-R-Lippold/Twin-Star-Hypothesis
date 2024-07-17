# Twin-Star-Hypothesis Repository
This repository contains the source code for the MATLAB-based computational work for the article "Twin-star hypothesis and cycle-free $d$-partitions of $K_{2d}$" by Matthew Fyfe, Steven Lippold, Mihai Staic, and Alin Stancu. More specifically, in the article, the theory had reduced a conjecture known as the Twin-Star Hypothesis down to the statement that every cycle-free edge 4-partition of the complete graph $K_8$ is weakly equivalent to the Twin-Star Graph with 8 vertices or an edge partition such that one of its graphs in the edge partition is the graph $T_19$. Using the work in this repository, it was shown computationally that these two cases are weakly equivalent, which proves the Twin-Star Hypothesis in the case of cycle-free edge 4-partitions of $K_8$.

## MATLAB Files
The MATLAB work was conducted in three steps, which is given in this repository as three .m files:
1. Partition Generation: This first step generates all of the relevant partitions. Keep in mind we are considering these partitions up to weak equivalence, given the Symmetric Group actions on edge 4-partitions.
2. Klein Group Action: This second step uses a group action by the Klein Group to reduce the relevant partitions, up to the group action.
3. Check for $TS_4$: This last step checks to see all of the partitions who are at most three involutions away from having the twin-star graph as one of its graphs in the edge 4-partition. It can be verified by the MATLAB code this is true for all of them, so the output file should be empty.

A more in-depth description is given by the [TwinStarDescription](../TwinStarDescription.pdf) file. The files here are free for use under the MIT License provided in this repository.

## Reference Work
A preprint version of the article "Twin-star hypothesis and cycle-free $d$-partitions of $K_{2d}$" by Fyfe et al can be found at ----------------.
