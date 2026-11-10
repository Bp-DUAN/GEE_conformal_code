% Solving surface diffuse/ mean curvature flow by Conformal map
% Crank-Nicolson sheme predict X_h^{m+1/2} by symmetreic Euelr scheme
% n_Eule: number of using Eulerian steps
% set n_Euler>Steps then it reduces to the Euler scheme
% theta: represents \omega in the paper
% Sep.2026,Beiping DUAN
clear
n_Euler=1;theta=0;
%%%
% T=0.09;Steps=943;
% tau0=1e-4;tau1=1e-6;
% tau_list=[tau0*ones(1,879),tau1*ones(1,120)];
%%%
T=0.68;Steps=680;
t_list=T*([0:1:Steps]/Steps).^2;
tau_list=t_list(2:end)-t_list(1:end-1);
%%%%%%%%% Sphere mesh for EXAMPLE 5.1
% load('../mesh/Sphere_h3.mat');
% Points=Sphere_h3.POS; % surface at t=0
% Tris=Sphere_h3.TRIANGLES(:,1:3);
%%%%%%%%% Dumbbell_06/07 mesh for EXAMPLE 5.2
% load('../mesh/Dumbbell_07.mat');
% Points=Dumbbell_07.POS;
% Tris=Dumbbell_07.TRIANGLES(:,1:3);
%%%%%%%%%% box 1-1-16 mesh for EXAMPLE 5.5
load('../mesh/Box_1_1_16.mat');
Points= Box_1_1_16.POS;
Tris=Box_1_1_16.TRIANGLES(:,1:3);
%%%%%% Torus mesh for EXAMPLE 6
% load('../mesh/Torus_65_3_ref2_s.mat');
% Points=Torus_65_3_ref2_s.POS; % surface at t=0
% Tris=Torus_65_3_ref2_s.TRIANGLES(:,1:3);
% R=1;
% Vn=[(sqrt(sum(Points.^2,2))-R).*Points(:,1),...
%     (sqrt(sum(Points.^2,2))-R).*Points(:,2),...
%     Points(:,3).*(Points(:,1).^2+Points(:,2).^2)];
%%%%%% genus-2 mesh for EXAMPLE 7 
% load('../mesh/genus2.mat');
% Points=genus2.POS; % surface at t=0
% Tris=genus2.TRIANGLES(:,1:3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lp=length(Points);
[~,S_phi_0,~,~,~,~] = Comp_Surface_MS_R(Points,Tris);
Num_U=zeros(3*Lp,Steps+1);
Num_U(:,1)=reshape(Points,[3*Lp,1]);
Energy=zeros(Steps+1,1);Vol=zeros(Steps+1,1);
Energy(1)=Comp_surface_area(Points,Tris);
Vol(1)=Comp_Volume(Points,Tris);
mesh_quality=zeros(Steps+1,1);
mesh_q_list = get_quality(reshape(Num_U(:,1),[Lp,3]),Tris);
mesh_quality(1)=max(mesh_q_list);
t=0;
for n=1:Steps
    tau=tau_list(n);
    t=t+tau;
    [cell_normal,~,~] = Comp_normal(reshape(Num_U(:,n),[Lp,3]),Tris);
    [M_phi,S_phi_n,~,~,~,~] = Comp_Surface_MS_R(reshape(Num_U(:,n),[Lp,3]),Tris);
    M_phi=diag(sum(M_phi,2));
    WS_phi=theta*S_phi_0+(1-theta)*S_phi_n;
%     MCF_SD=M_phi;% uncomment this line and comment the next line for mean curvature flow
    MCF_SD=S_phi_n;% surface diffusion
    M_phi_n=Comp_Mass_Cell_n(reshape(Num_U(:,n),[Lp,3]),Tris,cell_normal);
    if n<=n_Euler
        Total_F=[zeros(3*Lp,1);M_phi_n*Num_U(:,n)];
        u_H=[kron(eye(3),WS_phi),-M_phi_n';M_phi_n,MCF_SD*tau]\Total_F;
        Num_U(:,n+1)=u_H(1:3*Lp);
    else
        %%%%% prediction by symmetric Euler scheme
        tau_bf=(tau_list(n)+tau_list(n-1))/2;
        u_half_back=(Num_U(:,n)+Num_U(:,n-1))/2;
        Total_F=[zeros(3*Lp,1);M_phi_n*u_half_back];
        u_half_forward=[kron(eye(3),WS_phi),-M_phi_n';M_phi_n,MCF_SD*tau_bf]\Total_F;
        %%%%%% correction
        Num_U_Pre=u_half_forward(1:3*Lp);
        [cell_normal,~,~] = Comp_normal(reshape(Num_U_Pre,[Lp,3]),Tris);
        [M_phi,S_phi_nhalf,~,~,~,~] = Comp_Surface_MS_R(reshape(Num_U_Pre,[Lp,3]),Tris);
        M_phi=diag(sum(M_phi,2));
%         MCF_SD=M_phi;% mean curvature flow
        MCF_SD=S_phi_nhalf;% surface diffusion
        M_phi_n=Comp_Mass_Cell_n(reshape(Num_U_Pre,[Lp,3]),Tris,cell_normal);
        WS_phi=theta*S_phi_0+(1-theta)*S_phi_n;
        Total_F=[kron(-eye(3)/2,WS_phi)*Num_U(:,n);M_phi_n*Num_U(:,n)];
        %%%%%% solve linear equation
        u_H=[kron(eye(3)/2,WS_phi),-M_phi_n';M_phi_n,MCF_SD*tau]\Total_F;
        Num_U(:,n+1)=u_H(1:3*Lp);
    end
    %%%%%
    Energy(n+1)=Comp_surface_area(reshape(Num_U(:,n+1),[Lp,3]),Tris);
    Vol(n+1)=Comp_Volume(reshape(Num_U(:,n+1),[Lp,3]),Tris);
    mesh_q_list = get_quality(reshape(Num_U(:,n+1),[Lp,3]),Tris);
    mesh_quality(n+1)=max(mesh_q_list);
    %%%%%%%%%%%%%%%%%%%%%%
    fprintf('--Time=%8.6e\n',t)
    clf,patch('faces',Tris,'vertices',reshape(real(Num_U(:,n+1)),[Lp,3]),...
        'facecol','w','edgecol','k','LineWidth',0.05); 
    axis equal;
    view(3);cameratoolbar;drawnow
    axis off
end
%%
path='/Volumes/code/GF_WConformal_1/SD/CN_SEuler_th0';
% folder='h3_Steps_24';
% folder='dumbbell_07';
% folder='Torus_65_3_ref2_s_tau1e3';
folder='Box_1_1_16_tau1e3';
% folder='genus2_Steps_300';
file_name=[path,'/',folder];
mkdir(path,folder)

save([file_name,'/Num_U'],'Num_U','-v7.3')
save([file_name,'/mesh_quality'],'mesh_quality')
save([file_name,'/Energy'],'Energy')
save([file_name,'/Tris'],'Tris')
save([file_name,'/Vol'],'Vol')
save([file_name,'/tau_list'],'tau_list')