EXERCISE 6

Functions description : 
    1. calc_target_length
    
    2. split_long_edges
    
    3. collapse_short_edges
    
    4. equalize_valences
        In this function we didn't want not lose too much time by flipping an edge forth and back.
        So we simulate the flip by changing the valence of the vertices involved (+1 on new edge's endpoint and -1 to the old edge's endpoints).
        We then compared the sum of squared valence deviation before and after this simulation and check if it was better to flip or not.
    
    5. tangential_relaxation
        For this function we used the formula found at page 27 of the '08 Remeshing' PDF. To do so, we compute the laplace point and to not lose
        time computing the matrix/vertex multiplication we did the whole development on paper and just compute the final result.
        We take a lambda value of 1 (it is a constant value used to control the smoothing) and then we update the position by adding the update vector.