function [cell_normal, vertex_normal,area_t] = Comp_normal(Points,Tris,varargin)
% Compute the outward normal vector for a surface given by (Points, Tris)
% cell_normal -- Lt*3 the cell normal vector
% vertex_normal -- Lp*3 the vertex normal vector
% Nov,2023,Beiping
Lp=length(Points);
Lt=length(Tris);
V_ij=Points(Tris(:,2),:)-Points(Tris(:,1),:);
V_ik=Points(Tris(:,3),:)-Points(Tris(:,1),:);
outer=cross(V_ij,V_ik,2);
area_t=sqrt(sum(outer.^2,2));
cell_normal=outer./kron(ones(1,3),area_t);
area_t=area_t/2;
vertex_normal=0;
if nargin==3
    vertex_normal=zeros(Lp,3);
    for n=1:Lt
        vertex_normal(Tris(n,1),:)=vertex_normal(Tris(n,1),:)+area_t(n)*cell_normal(n,:);
        vertex_normal(Tris(n,2),:)=vertex_normal(Tris(n,2),:)+area_t(n)*cell_normal(n,:);
        vertex_normal(Tris(n,3),:)=vertex_normal(Tris(n,3),:)+area_t(n)*cell_normal(n,:);
    end
    vertex_normal=vertex_normal./(kron(ones(1,3),sqrt(sum(vertex_normal.^2,2))));
end

