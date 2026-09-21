# Godot-AVBD
## what is AVBD?
AVBD means Augmented Vertex Block Descent. It's a position based physics framework based on VBD,using Augmented Lagrange method to scaling the pentration force , which is very basic in position based physics.

It can handle a lot of rigidbodies cause its small size matrix of hessian and constraint between the huge and the small.


## what is the benefit
Of course it does , first , it can handle more objects and more accurate than those physics engine who use impulse projection method

And in the future it is possible to move all the compute to GPU and some company such as nvidia already did it.


## how does AVBD works?

### core formula
AVBD is based on Integration Method of Hessian based Newton Formula(Extended of course)
Then we can have 
1. the force formula which is H_i*delta_x_i = f_i
2. the Hessian Iteration formula which is H_i = M_i / (delta_t)^2 + Sum(H_ij)

transform the formula 1 to this:

3. delta_x_i = Inverse((H_i)) * f

### physics basic
as we know , we pursue minimum energy in realtime physics sim.
so we have this:

E_j(x) = 1/2 * k_j(C_j(x))^2

k_j is some kind of stiffness representative

C_j is the those kind of "penalty function" , which can tell you how close you are to the physical constraint

### avbd magic!
so usually here , xpbd use a simple lagrange multipliers to get them in position , avbd use lagrange multipliers,too.

but he use Augmented Lagrange which is the one who can achieve almost infinite stiffness.

so we have:

E_j(x) = 1/2* k_j(C_j(x))^2 + lambda_j * C_j(x)

lambda_j is the lagrange multipliers!

### so what's next
this is not those kind of tutorial article so that's all I wanna show you if you are a physics engine dev you know what I know

[if you really really want to know](v2-6e8c45cac4eb0cacfc6b952955f36cae_1440w.jpg)

