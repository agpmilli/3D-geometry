Exercise 2 :

	LaplacianSmoothing() :
		1) We instantiate two variables designed to be the length of the curve before (L) and after (L1) the transformation.
		2) We create a temporary matrix (points_tmp) containing the transformed points.
		3) We iterate through the points and compute the length between the current point and his next neighbor (add it to L).
		4) We compute the new point's coordinates using the formula given in the exercise and we store it in the temporary matrix (points_tmp)
		5) We compute the new point's coordinates for points with indexes 0 and num_points-1 because they are a bit special.
		6) We compute the new length of the curve by iterating through the temporary matrix. (add it to L1)
		7) At the end we replace the old points by the new ones and uniformly scale the curve back to its original length using L and L1. (we multiply the new points by L/L1)
		
	
	OsculatingCircle() :
		1) We instantiate two variables designed to be the length of the curve before (L) and after (L1) the transformation.
		2) We create a temporary matrix (points_tmp) containing the transformed points.
		3) We iterate through the points and compute the length between the current point and his next neighbor (add it to L).
		4) During this iteration we also compute the center of the circumscribed circle created by the 3 points with indexes i-1, i and i+1 (i being the current point).
		5) We computer the new point's coordinates using the formula given in the exercise and we store in the temporary matrix. (points_tmp)
		6) We compute the new length of the curve by iterating through the temporary matrix. (add it to L1)
		7) At the end we replace the old points by the new ones and uniformly scale the curve back to its original length using L and L1. (we multiply the new points by L/L1)
		
Exercise 3 :
	
	ReconstructCurveStep() :
		1) We compute a variable s by multiplying epsilon and n.
		2) We know that phi(s) is the integral of k(s) with respect to s. So we compute phi(s) for each k(s) given in the exercise.
		3) We compute the new points by adding to the precedent point the tangent.
		4) At the end we compute the tangent for next run using the formula (tangent = Vec2(cos(phi), sin(phi))).
		5) We comment/uncomment the wanted phi in our method and the wanted c(0) in the "MainWindow" to get it work.
		