% Solving solid-state dewetting model by Conformal map
% Crank-Nicolson sheme
% n_Eule: number of using Eulerian steps
% set n_Euler>Steps then it reduces to the Euler scheme
% theta: represents \omega in the paper
% Sep.2026,Beiping DUAN
clear
T=1;Steps=200; 
theta=1;
t_list=T*([0:1:Steps]/Steps);
tau_list=t_list(2:end)-t_list(1:end-1);
alpha=2/3*pi;beta=100;
n_Euler=1;
%%%%%%%%%%%% below is the initial mesh of EXAMPLE 5.3
load('./half_ellipsoid_levels/half_ellipsoid_l4.mat')
Tris_bound=load('./half_ellipsoid_levels/ellipsoid_4_boundary.txt');
Points=half_ellipsoid_l4.POS;
Tris=half_ellipsoid_l4.TRIANGLES(:,1:3);
%%%%%%%%%%%%% below is the initial mesh of EXAMPLE 5.4
% load('./Box_32_32_1_B_r.mat')
% Points=Box_32_32_1_B.POS;
% Tris=Box_32_32_1_B.TRIANGLES(:,1:3);
% Tris_bound=Box_32_32_1_B.LINES(:,1:2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Points_bound=Points(Tris_bound(:,1),:);
dim_bound=length(Points_bound);
Lp=length(Points);
Index=[1:Lp]';
Index(Tris_bound(:,1))=[];
[~,S_phi_0,~,~,~,~] = Comp_Surface_MS_R(Points,Tris);

Num_U=zeros(3*Lp,Steps+1);
Num_U(:,1)=reshape(Points,[3*Lp,1]);
Curve_area=zeros(Steps+1,1);
Surface_area=zeros(Steps+1,1);
Vol=zeros(Steps+1,1);

mesh_quality=zeros(Steps+1,1);
mesh_q_list = get_quality(reshape(Num_U(:,1),[Lp,3]),Tris);
mesh_quality(1)=max(mesh_q_list);
M_conv = spdiags(0.5*ones(dim_bound,1), 1, dim_bound, dim_bound) ...
       - spdiags(0.5*ones(dim_bound,1), -1, dim_bound, dim_bound);
M_conv(dim_bound,1)=0.5;M_conv(1,dim_bound)=-0.5;
[I1, J1, V1]=find(M_conv);           % 只取非零元素的行、列、值
B_term2=sparse(Tris_bound(I1,1), Tris_bound(J1,1), V1, Lp, Lp);
II=[0,1,0;-1,0,0;0,0,0];
t=0;
B_term1=sparse(Lp*3,Lp*3);
for n=1:Steps
    tau=tau_list(n);
    t=t+tau;
    u_H_Euler=zeros(4*Lp,1);
    u_half_forward=zeros(4*Lp,1);
    u_H_CN=zeros(4*Lp,1);
    %%%%%% boundary
    Points_n=reshape(Num_U(:,n),[Lp,3]);
    P_bound_n=Points_n(Tris_bound(:,1),1:2);
    %%%%%% Energy Volume
    Curve_area(n)=Comp_area(P_bound_n);
    Surface_area(n)=Comp_surface_area(Points_n,Tris);
    Vol(n)=Comp_Volume(Points_n,Tris);
    %%%%%%%%%%%%%%%%%%%%
    [tan_curve,~]=Comp_curve_tan(P_bound_n);
    cell_n_bound=[tan_curve(:,2),-tan_curve(:,1)];
    Mnxnx_curve=Comp_curve_M_weight(P_bound_n,cell_n_bound(:,1).^2);
    Mnxny_curve=Comp_curve_M_weight(P_bound_n,cell_n_bound(:,1).*cell_n_bound(:,2));
    Mnyny_curve=Comp_curve_M_weight(P_bound_n,cell_n_bound(:,2).^2);
    [I2, J2, V2]=find(Mnxnx_curve);           % 只取非零元素的行、列、值
    [I3, J3, V3]=find(Mnxny_curve);           % 只取非零元素的行、列、值
    [I4, J4, V4]=find(Mnyny_curve);           % 只取非零元素的行、列、值
    M_nxnx=sparse(Tris_bound(I2,1), Tris_bound(J2,1), V2, Lp, Lp);
    M_nxny=sparse(Tris_bound(I3,1), Tris_bound(J3,1), V3, Lp, Lp);
    M_nyny=sparse(Tris_bound(I4,1), Tris_bound(J4,1), V4, Lp, Lp);
    B_term1(1:2*Lp,1:2*Lp)=[M_nxnx,M_nxny;M_nxny',M_nyny];
    %%%%%% surface
    [cell_n,~,~] = Comp_normal(Points_n,Tris);
    [M_phi,S_phi_n,~,~,~,~] = Comp_Surface_MS_R(Points_n,Tris);
    M_phi=diag(sum(M_phi,2));
    MCF_SD=S_phi_n;
    M_phi_N=Comp_Mass_Cell_n(reshape(Num_U(:,n),[Lp,3]),Tris,cell_n);
    WS_phi=theta*S_phi_0+(1-theta)*S_phi_n;
    if n<=n_Euler
       Total_F=[(cos(alpha)/2*kron(II,B_term2)+...
                1/beta/tau*B_term1)*Num_U(:,n);M_phi_N*Num_U(:,n)];
       M_Total=[kron(eye(3),WS_phi)-cos(alpha)/2*kron(II,B_term2)+...
                1/beta/tau*B_term1,-M_phi_N';M_phi_N,MCF_SD*tau];
       M_Total_r=M_Total;Total_F_r=Total_F;
       M_Total_r(2*Lp+Tris_bound(:,1),:)=[];
       M_Total_r(:,2*Lp+Tris_bound(:,1))=[];
       Total_F_r(2*Lp+Tris_bound(:,1))=[];
       u_H_Euler([1:2*Lp,2*Lp+Index',3*Lp+1:4*Lp]')=M_Total_r\Total_F_r;
       Num_U(:,n+1)=u_H_Euler(1:3*Lp);
    else
     %%%%%%%%%%%%%%%%%%%%%%%%%%%  Euler-CN %%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%% prediction by SEuler scheme
       tau_bf=(tau_list(n)+tau_list(n-1))/2;
       u_half_back=(Num_U(:,n)+Num_U(:,n-1))/2;
       Total_F=[(cos(alpha)/2*kron(II,B_term2)+...
                1/beta/tau_bf*B_term1)*u_half_back;M_phi_N*u_half_back];
       M_Total=[kron(eye(3),WS_phi)-cos(alpha)/2*kron(II,B_term2)+...
                1/beta/tau_bf*B_term1,-M_phi_N';M_phi_N,MCF_SD*tau_bf];
       M_Total_r=M_Total;Total_F_r=Total_F;
       M_Total_r(2*Lp+Tris_bound(:,1),:)=[];
       M_Total_r(:,2*Lp+Tris_bound(:,1))=[];
       Total_F_r(2*Lp+Tris_bound(:,1))=[];
       u_half_forward([1:2*Lp,2*Lp+Index',3*Lp+1:4*Lp]')=M_Total_r\Total_F_r;
       Points_nhalf=reshape(u_half_forward(1:3*Lp),[Lp,3]);
       %%%%%%%%%%% correction by the CN scheme
       %%%%%%%%%% Boundary
       P_bound_nhalf=Points_nhalf(Tris_bound(:,1),1:2);
       [tan_curve,~]=Comp_curve_tan(P_bound_nhalf);
       cell_nhalf_bound=[tan_curve(:,2),-tan_curve(:,1)];
       Mnxnx_curve_half=Comp_curve_M_weight(P_bound_nhalf,cell_nhalf_bound(:,1).^2);
       Mnxny_curve_half=Comp_curve_M_weight(P_bound_nhalf,cell_nhalf_bound(:,1).*cell_nhalf_bound(:,2));
       Mnyny_curve_half=Comp_curve_M_weight(P_bound_nhalf,cell_nhalf_bound(:,2).^2);
       [I2, J2, V2]=find(Mnxnx_curve_half);           % 只取非零元素的行、列、值
       [I3, J3, V3]=find(Mnxny_curve_half);           % 只取非零元素的行、列、值
       [I4, J4, V4]=find(Mnyny_curve_half);           % 只取非零元素的行、列、值
       M_nxnx_half=sparse(Tris_bound(I2,1), Tris_bound(J2,1), V2, Lp, Lp);
       M_nxny_half=sparse(Tris_bound(I3,1), Tris_bound(J3,1), V3, Lp, Lp);
       M_nyny_half=sparse(Tris_bound(I4,1), Tris_bound(J4,1), V4, Lp, Lp);
       B_term1(1:2*Lp,1:2*Lp)=[M_nxnx_half,M_nxny_half;M_nxny_half',M_nyny_half];
    %%%%%%%%% surface
       [cell_nhalf,~,~] = Comp_normal(Points_nhalf,Tris);
       [M_phi,S_phi_nhalf,~,~,~,~] = Comp_Surface_MS_R(Points_nhalf,Tris);
       M_phi=diag(sum(M_phi,2));
       MCF_SD=S_phi_nhalf;
       WS_phi=theta*S_phi_0+(1-theta)*S_phi_n;
       M_phi_N=Comp_Mass_Cell_n(Points_nhalf,Tris,cell_nhalf);
       Total_F=[(kron(-eye(3)/2,WS_phi)+cos(alpha)/2*kron(II,B_term2)+...
                1/beta/tau*B_term1)*Num_U(:,n);M_phi_N*Num_U(:,n)];
       M_Total=[kron(eye(3)/2,WS_phi)-cos(alpha)/2*kron(II,B_term2)+...
         1/beta/tau*B_term1,-M_phi_N';M_phi_N,MCF_SD*tau];
     
       M_Total_r=M_Total;Total_F_r=Total_F;
       M_Total_r(2*Lp+Tris_bound(:,1),:)=[];
       M_Total_r(:,2*Lp+Tris_bound(:,1))=[];
       Total_F_r(2*Lp+Tris_bound(:,1))=[];
       u_H_CN([1:2*Lp,2*Lp+Index',3*Lp+1:4*Lp]')=M_Total_r\Total_F_r;
       Num_U(:,n+1)=u_H_CN(1:3*Lp);
    %%%%%
    end
    mesh_q_list = get_quality(reshape(Num_U(:,n+1),[Lp,3]),Tris);
    mesh_quality(n+1)=max(mesh_q_list);
    %%%%%%%%%%%%%%%%%%%%%%
    fprintf('--Time=%8.6e\n',t)
    clf,patch('faces',Tris,'vertices',reshape(real(Num_U(:,n+1)),[Lp,3]),...
        'facecol','w','edgecol','k','LineWidth',0.1); 
    axis equal;
%     axis off;
    view(3);cameratoolbar;drawnow
end
%%%%%% Energy Volume
Points_n=reshape(Num_U(:,n+1),[Lp,3]);
P_bound_n=Points_n(Tris_bound(:,1),1:2);
Curve_area(n+1)=Comp_area(P_bound_n);
Surface_area(n+1)=Comp_surface_area(Points_n,Tris);
Vol(n+1)=Comp_Volume(Points_n,Tris);
%%
path='/Volumes/code/GF_WConformal_1/SD_b/CN_SEuler_th1/';
folder='half_ellipsoid_l2_Steps_200_t1/alpha_23_beta_100';
% folder='half_Box_32_32_1_steps_200_t01/alpha23_beta_100';
file_name=[path,'/',folder];
mkdir(path,folder)

save([file_name,'/Num_U'],'Num_U','-v7.3')
save([file_name,'/mesh_quality'],'mesh_quality')
save([file_name,'/Curve_area'],'Curve_area')
save([file_name,'/Surface_area'],'Surface_area')
save([file_name,'/Vol'],'Vol')
save([file_name,'/Tris'],'Tris')
save([file_name,'/Tris_bound','Tris_bound'])
save([file_name,'/tau_list'],'tau_list')