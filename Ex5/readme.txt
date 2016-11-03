EXERCISE 5.1
For curves we saw that under curvature flow any curve goes to a convex shape and then
converges to a point.
Do you experience an analogous behavior for surfaces? Experiment with various meshes,
time steps, and number of iterations. Briefly comment on your observations (no formal
proof expected).

When using the implicit smoothing on a surface, the max-planck for instance, there are no changes to be observed,
 even when varying the time steps in iterations.

EXERCISE 5.2
Iterate your method on the three provided cylinders. One of them shows a behavior that
is different from the other two. Can you explain what is happening? Is the result consistent
with the goal of the minimal surface optimization?

When applying the algorithm, we minimize any surface area, that means that any changes will lead
to an increase of surface area. After enough iterations, hte area will be minimal on the whole surface¨
and hence will not change anymore under further iterations.

The main difference between the 3 objects is that compared to the cylinder, the two other objects converge to a final and minimized form 
whereas the cylinder idsappears.

See screenshots in "Screenshots" folder