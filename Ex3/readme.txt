Group members: Thomas Batschelet, Alain Milliet, Raphaël Steinmann
Each one of us did 1/3 of the work.

1: Theory Exercise

The answer for the theory exercise is written in the pdf file "ex3-theory.pdf"

2: Coding exercise

2.1: Vertex valence of a triangle mesh
* First we iterate over all vertices in the mesh using a "vertex iterator"
* For each vertex, we iterate over all its neighbors (using a "vertex around vertex iterator") and for each neighbor we increment the "valence" variable by one.
	- Don't forget to set the "valence" variable back to 0 before iterating of the neighbors of the next vertex.
* Store the valence in the "vertex_valence" array.

2.2: Computing Vertex Normals

Here instead of explicitly declaring and using iterators as in ex 2.1, we use "for each" loops (though it does exactly the same thing, we found it more comfortable to use)

Compute Normals with constant weights:
* iterate over all vertices
* for each vertex v:
	* initialize a sum variable to 0
	* iterate over all triangles T (Face objects) incident to v:
		* for each triangle, increment the sum by the weight (here, weight=1) times the normal vector of T (which we obtained with the compute_face_normal function).
	* normalize the sum and store it in v_cste_weights_n array.
	
Compute Normals by area weights:
* here the process is exactly the same as above, except for the weights. The weights are now The area of the current triangle T instead of constantly being 1.
* We added a third for loop which stores all three vertices of the triangle T in an array called "corners" and then we pass this array to the "computeArea" function described under.
* ComputeArea basically takes in input an array contaning the three vertices of a triangle and returns the area of the triangle. See code and comments for more details.
* Inside the function, we use an assert statement to make sure the shape is indeed a triangle (i.e. it has three vertices).

Compute Normals with Angle Weights:
* This function works exactly the same way as the one with area weights except that we compute the weights with another function called computeAngle.
* Given the three vertices of the triangle, this function returns the angle of the vertex being in the center of the neighborhood. See code/comments for more details.
* again we use two assertions: one as before to check that the shape is indeed a triangle (ie has 3 vetices) and another one to ensure that only one of the three vertices is the 
  vetex in the center of the neighborhood.
