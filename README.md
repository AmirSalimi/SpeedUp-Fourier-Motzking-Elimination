# SpeedUp-Fourier-Motzking-Elimination
1) Here A is the matrix as Halfspace representation
I put a sample A matrix

2) The variable you want to eliminate first, should be the first e elements

3) e should be specified as the number of variables you want to eliminiate

This algorithm, at each stage of FME, runs another LP to check the feasibility with new hyperplanes
If it a generated Hyperplane does not have feasible solution then it removes it

This kind of algorithm, amortize the cost since each extra hyperplane  (if is redundant) can cause exponensially number of different hyperplanes in next stages. Think about it!!!
