function [Mass_n] = Comp_Mass_Cell_n(Points,Tris,Cell_n,varargin)
% Compute the mass matrix \int \vu \cdot(\vn \xi)
% where n is the cell normal
% get mass-lumped version without setting varargin  
% Jan 2024 by Beiping
V_ij=Points(Tris(:,2),:)-Points(Tris(:,1),:);
V_ik=Points(Tris(:,3),:)-Points(Tris(:,1),:);
outer=cross(V_ij,V_ik,2);
area_t=sqrt(sum(outer.^2,2))/2;

Vx=area_t.*Cell_n(:,1);
Vy=area_t.*Cell_n(:,2);
Vz=area_t.*Cell_n(:,3);

S1=[Tris,Tris,Tris];
S2=kron(Tris,ones(1,3));

ijk=[2,1,1,1,2,1,1,1,2];
Mx=sparse(S1,S2,kron(Vx/12,ijk));
My=sparse(S1,S2,kron(Vy/12,ijk));
Mz=sparse(S1,S2,kron(Vz/12,ijk));

if nargin==3
    Lp=length(Points);
    Mx=spdiags(sum(Mx,2),0,Lp,Lp);
    My=spdiags(sum(My,2),0,Lp,Lp);
    Mz=spdiags(sum(Mz,2),0,Lp,Lp);
end
Mass_n=[Mx,My,Mz];
end

