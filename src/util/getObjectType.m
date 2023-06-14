function objectType = getObjectType(obj)
% getObjectType Returns the object type of input object.
% Inputs
%     obj     a initialized object of class SuperQuadrics or SuperEllipse
% Outputs
%     objectType  a string value of objecty type, including sphere,
%     ellip(ellipsoid or ellipse), superquadric

if size(obj.eps, 2)==1
%     object is 2D 
    dimension = 2;
else
%     object is 3D 
    dimension = 3;
end

if isequal(obj.eps, ones(1, dimension-1)) 
    objectType='ellip';
    if isequal(obj.a./obj.a(1), ones(1,dimension))
        objectType='sphere';
    end
else 
    objectType = 'superquadrics';
end

end