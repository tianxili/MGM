function adj_norm = edgenorm(fitlist)
% Function to calculate the group L2 norm for each pair of edges
% Input: 
% fitlist_post: the fitted parameter path
%         fitlist = struct('lambda_j',zeros(q,L),'lambda_jk',zeros(q1,L),'ita_0',zeros(p,L),'ita_j',zeros(p,q,L),...
%                          'phi_0', zeros(p,p,L), 'phi_j', zeros(p,p,q,L));
% Output: 
% adj_norm
[p p q L] = size(fitlist.phi_j);
adj_norm = struct('zz', zeros(q,q,L), 'zy', zeros(q,p,L), 'yy', zeros(p,p,L));
%% edges between z_j and z_k
for j = 1:q-1
    for k = j+1:q
        ind = transind(q,j,k);
        adj_norm.zz(j,k,:) = fitlist.lambda_jk(ind,:).^2;
        adj_norm.zz(k,j,:) = adj_norm.zz(j,k,:);
    end
end
%% edges between z_j and y_k
for j = 1:q
    for k = 1:p
        par = zeros(p+1,L);
        par(1,:) = fitlist.ita_j(k,j,:);
        par(2:p+1,:) = reshape(fitlist.phi_j(:,k,j,:),p,L);
        adj_norm.zy(j,k,:) = sum(par.^2,1);
    end
end
%% edges between y_j and y_k
for j = 1:p-1
    for k = j+1:p
        par = zeros(1+q, L);
        par(1,:) = fitlist.phi_0(j,k,:);
        par(2:q+1,:) = reshape(fitlist.phi_j(j,k,:,:),q,L);
        adj_norm.yy(j,k,:) = sum(par.^2,1);
        adj_norm.yy(k,j,:) = adj_norm.yy(j,k,:);
    end
end
clear par tmp;
end

        
        
        
            
        
