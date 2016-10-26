1.1 Uniform Laplace Curvature
TO do this we iterate over all vertices v, check if v is a boundary vertex. If that's the case we proceed to the next vertex. If v's not a boundary vertex, we iterate
over all its neighbors to calculate the sum of the uniform Laplace approximation formula. Once the sum is computed, we just have to calculate its norm and divide it by N.
N is the valence of v, which we computed before iterating over v's neighbors.
Then we just store the result in the v_unicurvature property array.

1.2 Laplace-Beltrami Curvature
- fill in the e_weight and v_weight property arrays with the calc_weights() function
- iterate over each vertex v of the mesh
- iterate over every neighbor vi of v, if vi is a boundary vertex proceed to the next neighbor
- get the edge e between v and vi with "find_edge" and the the edge-weight of e from the correct weight property array of the mesh
- compute the sum as expressed in the laplace beltrami operator formula. At the end we multiply the sum by w (as in the formula) which is the vertex-weight of v
- store the norm of the sum in the v_curvature property array

1.3
TODO

2.1 Uniform_smooth
We iterate over the non-boundary vertices of the mesh
For each vertex, we compute its valence and the sum of its neighbors position
We then compute the vector that will be added to its current position, being : (sum(neighbors) / valence) - current_position
We add this vector to the current position to find the new position

2.2 smooth
We precompute the edge weights using calc_edges_weights()
We iterate over the non-boundary vertices of the mesh
For each vertex, we iterate over its neighbors and compute the edge between the current point and the current neighbor and find the weight of this edge
We then compute the sum of these weights and the sum of this weights * vector between the neighbor and the current point.
We compute the vector that will be added to the current position, being : 1/sum(wi) * sum(wi * vec(neighbor, v))
Finally we had half of the result to the current position to get the new position of current point

3.1 uniform_laplacian_enhance_feature
First of all we iterate over the vertices and put the position of each vertex in a std::vector
We do the transformation on the vertices
Finally we reiterate over the vertices and adjust the position with the formula :
position = after + alpha * (before - after)
with after being the current position.

3.2 laplace_beltrami_enhance_feature
First of all we iterate over the vertices and put the position of each vertex in a std::vector
We do the transformation on the vertices
Finally we reiterate over the vertices and adjust the position with the formula :
position = after + alpha * (before - after)
with after being the current position.

Image 1 : bunny -> 2 x enhancement laplace_beltrami(n_iter = 10, alpha = 5) -> 1 x enhancement uniform_laplacian(n_iter=10, alpha=2)
-> 2 x smooth laplace_beltrami(n_iter=10)

Image 2 :

Image 3 :


From a signal processing point of view, what are the effects of the operations you apply, and why do they produce the results you show?