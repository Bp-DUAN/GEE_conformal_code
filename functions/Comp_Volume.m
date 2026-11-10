function [Volume] =Comp_Volume(Points,Tris)
% Compute the volume of a polyhedron
% one needs to use 'Tris_to_tris' first 
% to ensure the cell normal is outward
% Nov,2023,Beiping
Lt=length(Tris);
[cell_normal,~,area_t]=Comp_normal(Points,Tris);
Volume=0;
for n=1:Lt
    Volume=Volume+cell_normal(n,:)*(Points(Tris(n,1),:)...
        +Points(Tris(n,2),:)+Points(Tris(n,3),:))'*area_t(n);
end
Volume=Volume/9;
