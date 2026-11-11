function [mesh_q] = get_quality(Points,Tris)
% get the mesh quality of given mesh
% mesh_q=circumcircle/inscribed circle
% Lt=length(Tris);
% Oct. 2022 by Beiping 
[Lt,~]=size(Tris);
mesh_q=zeros(Lt,1);
for n=1:Lt
    p=Points(Tris(n,:)',:);
%     R=sqrt(sum((sum(p,1)/3-p(1,:)).^2));
    S=sqrt(det([p(2,:)-p(1,:);p(3,:)-p(1,:)]*[p(2,:)-p(1,:);p(3,:)-p(1,:)]')); % double area
    L=sum([p(2,:)-p(1,:);p(3,:)-p(2,:);p(3,:)-p(1,:)].^2,2);
    L=sqrt(L);
    r=S/sum(L);% radius of the maximal inscribed circle
    R=prod(L)/(2*S);
    mesh_q(n)=R/r;
end

