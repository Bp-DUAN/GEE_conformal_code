%% convergence test in TIME of CN scheme for solid deweting 
%%%% for CN scheme
hht={'half_ellipsoid_l4_Steps_24'...
    'half_ellipsoid_l4_Steps_32',...
    'half_ellipsoid_l4_Steps_40',...
    'half_ellipsoid_l4_Steps_48',...
    'half_ellipsoid_l4_Steps_56'};
%%%% for Euler scheme
% hht={'half_ellipsoid_l4_Steps_64',...
%     'half_ellipsoid_l4_Steps_128'...
%     'half_ellipsoid_l4_Steps_256',...
%     'half_ellipsoid_l4_Steps_512',...
%     'half_ellipsoid_l4_Steps_1024'};
stheta={'th1/','th05/','th0/'};
l_theta=[1,0.5,0];
Nt=length(hht);
Error_u_max=zeros(1,Nt-1);
Order_u_max=zeros(1,Nt-2);
Points_n_CN=cell(1,Nt);
Tris_CN=cell(1,Nt);
ERROR_d8=cell(length(l_theta),2);
for kth=1:length(l_theta)
    for k=1:Nt
        path_CN=['/Volumes/code/GF_WConformal_1/SD_b/CN_SEuler_',stheta{kth},hht{k},'_t02/alpha_12_beta_100'];
%         path_CN=['/Volumes/code/GF_WConformal_1/SD_b/Euler_',stheta{kth},hht{k},'_t02/alpha_23_beta_100'];
        folder_CN='/Num_U';
        folder_CN_t='/Tris';
        load([path_CN,folder_CN])
        load([path_CN,folder_CN_t])
        Lp=size(Num_U,1)/3;
        Points_n_CN{k}=reshape(Num_U(:,end),[Lp,3]);
        Tris_CN{k}=Tris;
    end
    for k=1:Nt-1
        P=Points_n_CN{k};
        Q=Points_n_CN{k+1};
        %%%%%%%%
        [dist_p2t, proj_p2t] = point2trimesh( ...
        'Faces',        Tris, ...
        'Vertices',     Q, ...
        'QueryPoints',  P);
        [dist_q2t, proj_q2t] = point2trimesh( ...
        'Faces',        Tris, ...
        'Vertices',     P, ...
        'QueryPoints',  Q);
        Error_u_max(k)=(max(abs(dist_p2t))+max(abs(dist_q2t)))/2;
        fprintf('--theta=%3.1f, k=%d\n',l_theta(kth),k)
    end
    for k=1:Nt-2
        %%%% for CN scheme
        Order_u_max(k)=log(Error_u_max(k)/Error_u_max(k+1))/log((8*(k+3))/(8*(k+2)));
        %%%% for Euler scheme
%         Order_u_max(k)=log(Error_u_max(k)/Error_u_max(k+1))/log(2);
    end
    ERROR_d8{kth,1}=Error_u_max;
    ERROR_d8{kth,2}=Order_u_max;
end
save('ERROR_d8','ERROR_d8')