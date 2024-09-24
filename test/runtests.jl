using BubbleSpaces
using Gridap
using Gridap.ReferenceFEs
using Test


N = 10;
m = CartesianDiscreteModel((0, 1, 0, 1), (N, N)) |> simplexify
labels = get_face_labeling(m)
add_tag_from_tags!(labels,"diri1",[6,])
add_tag_from_tags!(labels,"diri0",[1,2,3,4,5,7,8])

trian = Triangulation(m)
dx = Measure(trian, 3)
u0 = VectorValue(0,0)
u1 = VectorValue(1,0)

B_test = TestFESpace(m, BubbleRefFE(VectorValue{2, Float64}, TRI))
B_trial = TrialFESpace(B_test)
U_test = TestFESpace(m, LagrangianRefFE(VectorValue{2, Float64}, TRI, 1); conformity= :H1, dirichlet_tags= ["diri0","diri1"])
U_trial = TrialFESpace(U_test, [u0,u1])
Q = TestFESpace(m, LagrangianRefFE(Float64, TRI, 1); conformity= :H1, constraint= :zeromean)
P = TrialFESpace(Q)

Xb = MultiFieldFESpace([U_trial, B_trial, P])
Yb = MultiFieldFESpace([U_test, B_test, Q])

@inline a(u, v) = ∫(∇(u)⊙∇(v))dx 
@inline b(u, p) = ∫(∇⋅u * p)dx

@inline A((u, p), (v, q)) = a(u, v) + b(v, p) + b(u, q)
@inline AB((uc, ub, p), (vc, vb, q)) = a(uc, vc) + a(ub, vb) + a(uc, vb) + a(ub, vc) + b(uc, q) + b(ub, q) + b(vc, p) + b(vb, p)
@inline AB_2((uc, ub, p), (vc, vb, q)) = A((uc+ub, p), (vc+vb, q))
@inline l((v, q)) = 0.

op_b = AffineFEOperator(AB_2, l, Xb, Yb)
uch, ubh, ph = solve(LUSolver(), op_b)

A1 = assemble_matrix(AB, Xb, Yb);
A2 = assemble_matrix(AB_2, Xb, Yb);
@test A1 == A2
