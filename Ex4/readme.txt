1.1
TODO

1.2
TODO

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