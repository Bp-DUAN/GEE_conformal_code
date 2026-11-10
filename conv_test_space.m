%% convergence in SPACE   for mean curvature flow   %%%%%
Nx=4;
Error_u=zeros(1,Nx);
Order_u=zeros(1,Nx-1);
u=cell(1,Nx);
hh={'Sphere_h0','Sphere_h1','Sphere_h2','Sphere_h3'};
h=zeros(1,Nx);
for k=1:Nx
    path_CN=['/Volumes/code/GF_WConformal_1/MCF/CN_SEuler/',hh{k}];
    folder_CN='/Num_U';
    folder_CN_t='/Tris';
    load([path_CN,folder_CN])
    load([path_CN,folder_CN_t])
    Lp=size(Num_U,1)/3;
    Points_n_CN=reshape(Num_U(:,end),[Lp,3]);
    [max_h,L]=Get_length_Tris(reshape(Num_U(:,1),[Lp,3]),Tris);
    h(k)=max_h;
    Error_u(k)=max(abs(sqrt(sum(Points_n_CN.^2,2))-0.2*ones(Lp,1)));
end
for k=1:Nx-1
      Order_u(k)=log(Error_u(k)/Error_u(k+1))/log(h(k)/h(k+1));
end
%% convergence test in SPACE for solid deweting
Nx=5;
Error_u_max=zeros(1,Nx-1);
Order_u_max=zeros(1,Nx-2);
Points_n_CN=cell(1,Nx);
Tris_CN=cell(1,Nx);
hh={'half_ellipsoid_l0_Steps_200','half_ellipsoid_l1_Steps_200',...
    'half_ellipsoid_l2_Steps_200',...
    'half_ellipsoid_l3_Steps_200','half_ellipsoid_l4_Steps_200'};
h=zeros(1,Nx);
for k=1:Nx
    path_CN=['/Volumes/code/GF_WConformal_1/SD_b/CN_SEuler_th0/',hh{k},'_t02','/alpha_23_beta_100'];
    folder_CN='/Num_U';
    folder_CN_t='/Tris';
    load([path_CN,folder_CN])
    load([path_CN,folder_CN_t])
    Lp=size(Num_U,1)/3;
    Points_n_CN{k}=reshape(Num_U(:,end),[Lp,3]);
    Tris_CN{k}=Tris;
    [max_h,L]=Get_length_Tris(reshape(Num_U(:,1),[Lp,3]),Tris);
    h(k)=max_h;
end
for k=1:Nx-1
    P=Points_n_CN{k};
    Q=Points_n_CN{k+1};
    [dist_p2t, proj_p2t] = point2trimesh( ...
        'Faces',        Tris_CN{k+1}, ...
        'Vertices',     Q, ...
        'QueryPoints',  P);
    [dist_q2t, proj_q2t] = point2trimesh( ...
    'Faces',        Tris_CN{k}, ...
    'Vertices',     P, ...
    'QueryPoints',  Q);
    Error_u_max(k)=(max(abs(dist_p2t))+max(abs(dist_q2t)))/2;
end
for k=1:Nx-2
      Order_u_max(k)=log(Error_u_max(k)/Error_u_max(k+1))/log(h(k)/h(k+1));
end
save('Order_u_max','Order_u_max')