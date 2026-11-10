function [area] =Comp_surface_area(Points,Tris)
% Compute the surface area of a polyhedral surface in R^3
% Nov,2023, Beiping
V_ij=Points(Tris(:,2),:)-Points(Tris(:,1),:);
V_ik=Points(Tris(:,3),:)-Points(Tris(:,1),:);
outer=cross(V_ij,V_ik,2);
area=sum(sqrt(sum(outer.^2,2)))/2;

end

