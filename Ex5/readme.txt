EXERCISE 5.1
For curves we saw that under curvature flow any curve goes to a convex shape and then
converges to a point.
Do you experience an analogous behavior for surfaces? Experiment with various meshes,
time steps, and number of iterations. Briefly comment on your observations (no formal
proof expected).

When using the implicit smoothing on a surface, the max-planck for instance, there are no changes to be observed,
even when varying the time steps in iterations.This is because it isn't a boundary mesh, whereas the bunny
will see is form becoming convex until merging in a single point. If we augment the time steps, it will accelerate the process.

EXERCISE 5.2
Iterate your method on the three provided cylinders. One of them shows a behavior that
is different from the other two. Can you explain what is happening? Is the result consistent
with the goal of the minimal surface optimization?

When applying the algorithm, we minimize any local surface area, that means that any changes afterwards will lead
to an increase of local surface area. Another way to describe it, is that the mean curvature is equal to zero after 
every iteration. Since we have constraints on the boundaries vertices it will need to compensate the other normal vectors
to reach a zero mean value.

The main difference between the 3 objects is that compared to the cylinder, the two other objects converge to a final and minimized form 
whereas the cylinder will be able to shrink and keeping its zero mean curvature until it converges to a single point.

See screenshots in "Screenshots" folder