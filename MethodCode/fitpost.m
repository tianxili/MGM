function fitlist_post = fitpost(fitlist, method)
% Function to post process to fitted parameters by taking the maximum or
% minimum absolute values with signs. 
% Input:
% fitlist: an object from sepreg, including the solution path
% method: 'max' or 'min'
p = size(fitlist.ita_0,1);
[q L] = size(fitlist.lambda_j);
q1 = q*(q-1)/2;
fitlist_post = struct('lambda_j',zeros(q,L),'lambda_jk',zeros(q1,L),'ita_0',zeros(p,L),'ita_j',zeros(p,q,L),...
                      'phi_0', zeros(p,p,L), 'phi_j', zeros(p,p,q,L));
fitlist_post.lambda_j = fitlist.lambda_j;
fitlist_post.ita_0 = fitlist.ita_0;
if strcmp(method, 'max')
    [a ind] = max(abs(fitlist.lambda_jk),[],2);
    s = size(fitlist.lambda_jk);
    ind1 = sub2ind(s,repmat(1:s(1),1,L)',ind(:),reshape(repmat(1:L,s(1),1),L*s(1),1));
    fitlist_post.lambda_jk = reshape(fitlist.lambda_jk(ind1),q1,L);
    
    [a ind] = max(abs(fitlist.ita_j),[],3);
    s = size(fitlist.ita_j);
    ind1 = sub2ind(s,repmat(1:s(1),1,L*s(2))',repmat(reshape(repmat(1:s(2),s(1),1),s(1)*s(2),1),L,1),...
                    ind(:),reshape(repmat(1:L,s(1)*s(2),1),L*s(1)*s(2),1));
    fitlist_post.ita_j = reshape(fitlist.ita_j(ind1),p,q,L);
    
    [a ind] = max(abs(fitlist.phi_0),[],3);
    s = size(fitlist.phi_0);
    ind1 = sub2ind(s,repmat(1:s(1),1,L*s(2))',repmat(reshape(repmat(1:s(2),s(1),1),s(1)*s(2),1),L,1),...
                    ind(:),reshape(repmat(1:L,s(1)*s(2),1),L*s(1)*s(2),1));
    fitlist_post.phi_0= reshape(fitlist.phi_0(ind1),p,p,L);
    
    [a ind] = max(abs(fitlist.phi_j),[],4);
    s = size(fitlist.phi_j);
    ind1 = repmat(1:s(1),1,s(2)*s(3)*L)';
    ind2 = repmat(reshape(repmat(1:s(2),s(1),1),s(1)*s(2),1),s(3)*L,1);
    ind3 = repmat(reshape(repmat(1:s(3),s(1)*s(2),1),s(1)*s(2)*s(3),1),L,1);
    ind5 = reshape(repmat(1:L,s(1)*s(2)*s(3),1),s(1)*s(2)*s(3)*L,1);
    ind = sub2ind(s,ind1, ind2, ind3, ind(:),ind5);    
    fitlist_post.phi_j = reshape(fitlist.phi_j(ind),p,p,q,L);

elseif strcmp(method,'min')
    [a ind] = min(abs(fitlist.lambda_jk),[],2);
    s = size(fitlist.lambda_jk);
    ind1 = sub2ind(s,repmat(1:s(1),1,L)',ind(:),reshape(repmat(1:L,s(1),1),L*s(1),1));
    fitlist_post.lambda_jk = reshape(fitlist.lambda_jk(ind1),q1,L);
    
    [a ind] = min(abs(fitlist.ita_j),[],3);
    s = size(fitlist.ita_j);
    ind1 = sub2ind(s,repmat(1:s(1),1,L*s(2))',repmat(reshape(repmat(1:s(2),s(1),1),s(1)*s(2),1),L,1),...
                    ind(:),reshape(repmat(1:L,s(1)*s(2),1),L*s(1)*s(2),1));
    fitlist_post.ita_j = reshape(fitlist.ita_j(ind1),p,q,L);
    
    [a ind] = min(abs(fitlist.phi_0),[],3);
    s = size(fitlist.phi_0);
    ind1 = sub2ind(s,repmat(1:s(1),1,L*s(2))',repmat(reshape(repmat(1:s(2),s(1),1),s(1)*s(2),1),L,1),...
                    ind(:),reshape(repmat(1:L,s(1)*s(2),1),L*s(1)*s(2),1));
    fitlist_post.phi_0= reshape(fitlist.phi_0(ind1),p,p,L);
    
    [a ind] = min(abs(fitlist.phi_j),[],4);
    s = size(fitlist.phi_j);
    ind1 = repmat(1:s(1),1,s(2)*s(3)*L)';
    ind2 = repmat(reshape(repmat(1:s(2),s(1),1),s(1)*s(2),1),s(3)*L,1);
    ind3 = repmat(reshape(repmat(1:s(3),s(1)*s(2),1),s(1)*s(2)*s(3),1),L,1);
    ind5 = reshape(repmat(1:L,s(1)*s(2)*s(3),1),s(1)*s(2)*s(3)*L,1);
    ind = sub2ind(s,ind1, ind2, ind3, ind(:),ind5);    
    fitlist_post.phi_j = reshape(fitlist.phi_j(ind),p,p,q,L);
   
end
end

    
    
    
    


