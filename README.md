# Game-Theoretic Perspective on a Demand-Side Management Problem

Dear reader,  
I am thrilled to briefly guide you through this project, which deals with a game-theoretic perspective on a Demand-Side Management problem. Specifically, a large number of customers aim at satisfying their energy demand while minimizing their own procurement cost. The problem set takes into account various constraints and features high computational cost. In the [report](DSM_Game.pdf), you will find:  
- the conditions for which the general problem can be reduced to a convex-game framework;  
- an energetic interpretation of the problem via Lagrange Mechanics and the subsequent calculation of an overall potential function;  
- the computation of the unique Nash Equilibrium strategy profile via Projected Gradient Descent.  
- brief comments on the result  

As concerns the computation, the [`.m` file](MATLAB_Code_DSM.m) contains a short extract of the MatLab code where the Projected Gradient Descent is implemented. The other [`.csv` files](CSV_Archive.zip) contain the data of a specific instance of the problem (beware: they are necessary for the `.m` file to run properly).  

Thank you very much for your attention and I hope you will enjoy this project.  

Davide!


## Projected Gradient Descent Update
+![Linkedin Post GTC - Iteration Update.gif]()

## Nash Equilibrium Approaching vs Iteration Count
+![Linkedin Post GTC - Distance WRT NE.gif]()
