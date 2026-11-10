function [M_phi,S_phi,kappa,Area,min_tri,max_tri] = Comp_Surface_MS_R(Points,Tris)
% ------------------------------------------------------------
% To compute mass, stiffness matrix and curvature of a surface
% when use piecewise linear basis functions
% revised version, using sparse storage. 2022-5-28,Beiping
% ------------------------------------------------------------
Area=0;min_tri=1;max_tri=0;
Lt=length(Tris);
Lp=length(Points);
M1_vec=zeros(6*Lt,1);M2_vec=zeros(6*Lt,1);
Val_M=zeros(6*Lt,1);Val_S=zeros(6*Lt,1);
    for m=1:Lt
        i=Tris(m,1);j=Tris(m,2);k=Tris(m,3);
        ijk=sort([i,j,k]);
        i=ijk(1);j=ijk(2);k=ijk(3);
        M1_vec(6*(m-1)+1)=i;M2_vec(6*(m-1)+1)=j;
        M1_vec(6*(m-1)+2)=j;M2_vec(6*(m-1)+2)=k;
        M1_vec(6*(m-1)+3)=i;M2_vec(6*(m-1)+3)=k;
        M1_vec(6*(m-1)+4)=i;M2_vec(6*(m-1)+4)=i;
        M1_vec(6*(m-1)+5)=j;M2_vec(6*(m-1)+5)=j;
        M1_vec(6*(m-1)+6)=k;M2_vec(6*(m-1)+6)=k;      
        GXiU=[Points(j,:)-Points(i,:);Points(k,:)-Points(i,:)]';
        Metric_g=GXiU'*GXiU;
%         det_g=det(Metric_g);
        det_g=det(Metric_g);% to avoid arising coplex number in extreme cases
        sqrt_det_g=sqrt(det_g);Area=Area+sqrt_det_g;
        min_tri=min(min_tri,sqrt_det_g/2);
        max_tri=max(max_tri,sqrt_det_g/2);
        % Mass matrix
        Val_M(6*(m-1)+1)=1/24*sqrt_det_g;
        Val_M(6*(m-1)+2)=1/24*sqrt_det_g;
        Val_M(6*(m-1)+3)=1/24*sqrt_det_g;
        Val_M(6*(m-1)+4)=1/24*sqrt_det_g;
        Val_M(6*(m-1)+5)=1/24*sqrt_det_g;
        Val_M(6*(m-1)+6)=1/24*sqrt_det_g;
        % Stiffness Matrix
        Inv_g=[Metric_g(2,2),-Metric_g(1,2);...
            -Metric_g(2,1),Metric_g(1,1)]/det_g;
        Val_S(6*(m-1)+1)=[-1,-1]*Inv_g*[1;0]*sqrt_det_g/2;
        Val_S(6*(m-1)+2)=[1,0]*Inv_g*[0;1]*sqrt_det_g/2;
        Val_S(6*(m-1)+3)=[-1,-1]*Inv_g*[0;1]*sqrt_det_g/2;
        Val_S(6*(m-1)+4)=[-1,-1]*Inv_g*[-1;-1]*sqrt_det_g/4;
        Val_S(6*(m-1)+5)=[1,0]*Inv_g*[1;0]*sqrt_det_g/4;
        Val_S(6*(m-1)+6)=[0,1]*Inv_g*[0;1]*sqrt_det_g/4;
    end
    Area=Area/2;
    M_phi=sparse(M1_vec,M2_vec,Val_M);
    M_phi=M_phi+M_phi';
    S_phi=sparse(M1_vec,M2_vec,Val_S);
    S_phi=S_phi+S_phi';    
%     kappa=M_phi\(S_phi*Points);
    kappa=sparse(1:Lp,1:Lp,sum(M_phi,2))\(S_phi*Points);
end

