function [ v ] = unit_vector(V)

nV=norm(V);

v=zeros(size(V));
if nV>10^-4
    v = V/norm(V);
end

end

