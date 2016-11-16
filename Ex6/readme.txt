EXERCISE 6

Functions description : 
    1. calc_target_length
        - Average length :
		This part is quite trivial. We first iterate over every edge of the mesh to compute the mean edge length in the mesh. Then we iterate over every vertex v
		and set v's target_length property to mean_edge_length.
        
        - Height based :
        For the height based remeshing, we first of all get the minimum height of a vertex in the mesh. We then set the target_length of each vertex 
		to its height + the minimum height found before (we basically shift the domain to avoid negative target_length values).
		Finally, we scale the target_length of each vertex so that its mean value is equal to the user specified target length. We defined the user spec. target length
		as a local variable and gave it a value such that our remeshed Max Planck face looks approximately like the one in the pdf.
		Note: to get the height value of a vertex v we just use mesh_.position(v)[1].
        
        - Curvature-based adaptive remeshing :
        For adaptive remeshing to know what was the better curvature to use, we took a look at 
            
        Additionnal comments : 
        With our implementation we don't get the exact same result shown in your PDF but we get your result when we do the exact same thing without using the curvature 
        to compute the initial target length. It is pretty strange but it looks like the image you gave us for 'Adaptive Remeshing' comes when every target_length = -nan(ind).
        
    2. split_long_edges
    This method aims to split edges that are too long. We perform at most 100 iterations (we use a boolean variable to check if there's no more splitting to do
    and break if it's the case). At each iteration we check every edge e, compute it's target length from the target lengths of its two endpoints v0 and v1 and 
    compare it to its length to decide if e must be split or not. If the answer is yes, we add a new vertex v between v0 and v1 and split e in two edges. Then, 
    we compute v's normal and interpolate its target_length from the target lengths of v0 and v1, doing a simple average.
    
    3. collapse_short_edges
    This function will iterate over all edges and compare their length to the mean target length of their respective end vertices. If it is smaller than a certain threshold
    (i.e. 4/5 of the target length) it will merge the lower valence vertex into the higher one. We make sure, that the end vertices are either both or aren't boundary vertices
    and that the halfedges are collapsible in order to use the collapse function. If only one of them is boundary, we collapse the non boundary into the boundary vertex after
    checking if the halfedge is collapsible. Furthermore we turn the variable finished to false every time that a collapse occurs so that if nothing changes, the iteration 
    will stop and won't lead to unnecessary loops.
    
    4. equalize_valences
    In this function we didn't want to lose too much time by flipping an edge forth and back.
    So we simulate the flip by changing the valence of the vertices involved (+1 on new edge's endpoint and -1 to the old edge's endpoints).
    We then compared the sum of squared valence deviation before and after this simulation and check if it was better to flip or not.
    
    5. tangential_relaxation
    For this function we used the formula found at page 27 of the '08 Remeshing' PDF. To do so, we compute the Laplace point and to not lose
    time computing the matrix/vertex multiplication we did the whole development on paper and just compute the final result.
    We take a lambda value of 1 (it is a constant value used to control the smoothing) and then we update the position by adding the update vector.

