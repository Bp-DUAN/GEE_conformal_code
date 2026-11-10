function [max_L,L] = Get_length_Tris(Points,Tris)
% Return the maximal length of bars
% Nov. 2023, Beiping
    bars=[Tris(:,[1,2]);Tris(:,[1,3]);Tris(:,[2,3])];% Interior bars duplicated
    bars=unique(sort(bars,2),'rows');                % Bars as node pairs
    barvec=Points(bars(:,1),:)-Points(bars(:,2),:);  % List of bar vectors
    L=sqrt(sum(barvec.^2,2));                        % L = Bar lengths
    max_L=max(L);
end

