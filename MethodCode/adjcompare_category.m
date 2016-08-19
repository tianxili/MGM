function [sens spec total_p total_n sens_zz spec_zz total_p_zz total_n_zz sens_yy spec_yy total_p_yy total_n_yy sens_zy spec_zy total_p_zy total_n_zy ] = adjcompare_category(adj_thres, adj,p,q)
% Function to compare the estimated edgeset and the true edge set
% INPUT:
% adj_thres: 4 dim array (p+q)*(p+q)*L*H
% adj: the true adjascent matrix
% OUTPUT:
% sens: sensitivity L*H
% spec: specificity L*H
D = size(adj,1);
L = size(adj_thres, 3);
H = size(adj_thres, 4);
sens = zeros(L,H);
spec = zeros(L,H);
sens_zz = zeros(L,H);
spec_zz = zeros(L,H);
sens_yy = zeros(L,H);
spec_yy = zeros(L,H);
sens_zy = zeros(L,H);
spec_zy = zeros(L,H);
tmp = adj(triu(true(size(adj)),1));
zz_adj = adj(1:q,1:q);
tmp_zz = zz_adj(triu(true(size(adj(1:q,1:q))),1));
yy_adj = adj((q+1):(q+p),(q+1):(q+p));
tmp_yy = yy_adj(triu(true(size(adj((q+1):(q+p),(q+1):(q+p)))),1));
tmp_zy = adj(1:q,(q+1):(q+p));

total_p = sum(tmp);
total_n = (D)*(D-1)/2-total_p;

total_p_zz = sum(tmp_zz);
total_n_zz = (q)*(q-1)/2-total_p_zz;

total_p_yy = sum(tmp_yy);
total_n_yy = (p)*(p-1)/2-total_p_yy;

total_p_zy = sum(sum(tmp_zy));
total_n_zy = p*q-total_p_zy;


tmp = adj(triu(true(size(adj)),1));
for h = 1:H
    for l = 1:L
        tmp1 = adj_thres(:,:,l,h);
        tmp1 = tmp1(triu(true(size(tmp1)),1));
        tp = sum(tmp1.*tmp);
        tn = sum(tmp1+tmp==0);
        sens(l,h) = tp/total_p;
        spec(l,h) = tn/total_n;
        %% zz
        tmp1_zz = adj_thres(1:q,1:q,l,h);
        tmp1_zz = tmp1_zz(triu(true(size(tmp1_zz)),1));
        tp_zz = sum(tmp1_zz.*tmp_zz);
        tn_zz = sum(tmp1_zz+tmp_zz==0);
        sens_zz(l,h) = tp_zz/total_p_zz;
        spec_zz(l,h) = tn_zz/total_n_zz;
        %% yy
        tmp1_yy = adj_thres((q+1):(q+p),(q+1):(q+p),l,h);
        tmp1_yy = tmp1_yy(triu(true(size(tmp1_yy)),1));
        tp_yy = sum(tmp1_yy.*tmp_yy);
        tn_yy = sum(tmp1_yy+tmp_yy==0);
        sens_yy(l,h) = tp_yy/total_p_yy;
        spec_yy(l,h) = tn_yy/total_n_yy;
        %% zy
        tmp1_zy = adj_thres(1:q,(q+1):(q+p),l,h);
        tp_zy = sum(sum(tmp1_zy.*tmp_zy));
        tn_zy = sum(sum(tmp1_zy+tmp_zy==0));
        sens_zy(l,h) = tp_zy/total_p_zy;
        spec_zy(l,h) = tn_zy/total_n_zy;
    end
end



end

        
        
        
        
