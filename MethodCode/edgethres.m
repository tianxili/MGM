
function adj_thres = edgethres(adj_norm, thres)
% Function to obtain the adjascent matrix by thresholding the adj norm
% matrix

% Input:
% adjnorm: a structure with adj norm matrix zz zy yy
% K: length of thresholding vector
% Output: 
% adj_thres: a 4-dimentional array to record the adj matrix (p+q)*(p+q)*L*K
K = length(thres);
[q p L]= size(adj_norm.zy);
adj_thres = zeros(p+q, p+q, L, K);

for k = 1:K
    adj_thres(1:q, 1:q, :,k) = (adj_norm.zz > thres(k));
    adj_thres(q+1:end, q+1:end, :, k) = (adj_norm.yy >thres(k));
    adj_thres(1:q, q+1:end,:,k) = (adj_norm.zy>thres(k));
    for l = 1:L
    adj_thres(q+1:end, 1:q,l,k) = adj_thres(1:q, q+1:end,l,k)';
    end
end
end


    