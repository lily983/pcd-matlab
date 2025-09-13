function x_new = method6(obj, k)
%idea: 1, 2, 3, 4
% obj: superquadric
% k: scaling constant of new obj
% x_new: new encapsulating surface points 3 x m

obj_scale = SuperQuadrics({obj.a, obj,eps', [0,0], obj.tc, obj.q, obj.N});

x_new = obj_scale.GetPoints();

end