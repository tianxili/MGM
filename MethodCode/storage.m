
function fitlist = storage(fitlist, a0, beta, j, sigma2)
% Function to store the fitted parameter from one regression into the current parlist. 
% Input: 
% fitlist: current fitted parameter list with different copies
%           struct('lambda_j', zeros(q,L), 'lambda_jk', zeros(q1,2,L), 'ita_0', zeros(p,L), 'ita_j', zeros(p,q,2,L),... 
%                 'ita_jk', zeros(p,q1,3,L), 'phi_0', zeros(p,p,2,L), 'phi_j', zeros(p,p,q,3,L), 'phi_jk', zeros(p,p,q1,4,L));
% a0: fitted intercept L*1
% beta: fitted coefficient path dim*L
% type: current regression type 'z' (binary) or 'y' (continuous)
% j: current regression index
% sigma2: if exists, a scalor from LS regression to rescale the parameters
% Output:
% fitlist: restored fitlist by plugging in the current estimate (a0, beta).

p = size(fitlist.ita_0,1);
q = size(fitlist.lambda_j,1);
p1 = (p-1)*p/2;
L = size(beta,2); % length of path, same as lenght of tuning pars
start = 0;
if nargin ==4
    if(j>1)
    ind1 = transind(q,j,1:j-1);
    end
    if(j<q)
    ind2 = transind(q,j,j+1:q);
    end
    fitlist.lambda_j(j,:) = a0;
    if (j>1)
        fitlist.lambda_jk(ind1,2,:) = beta(start+1:start+(j-1),:);
    end
    start = start +j-1;
    if(j<q)
        fitlist.lambda_jk(ind2,1,:) = beta(start+1:start+(q-j),:);
    end
    start = start + q-j;
    fitlist.ita_j(:,j,1,:) = beta(start+1:start+p,:);
    start = start+p;
    %tmp = tril(ones(p)); %% if update phi_diagonal then no -1
    tmp = tril(ones(p),-1);
    tmp = reshape(repmat(tmp,1,L),p,p,L);
    tmp(tmp==1) = beta(start+1:start+p1,:);
    for l = 1:L
    fitlist.phi_j(:,:,j,1,l) = tmp(:,:,l)+tmp(:,:,l)'-diag(diag(tmp(:,:,l)));
    end
elseif nargin==5
    fitlist.ita_0(j,:)= a0./sigma2;
    beta = beta./repmat(reshape(sigma2,1,L),size(beta,1),1);
    fitlist.phi_0(j,j,1,:) = 1./sigma2;
    fitlist.phi_0(j,j,2,:) = 1./sigma2;
    fitlist.ita_j(j,:,2,:) = beta(start+1:start+q,:);
    start = start+q;
    if(j>1)
        fitlist.phi_0(j,1:j-1,2,:) = beta(start+1:start+j-1,:);
        fitlist.phi_0(1:j-1,j,2,:) = beta(start+1:start+j-1,:);
    end
    start = start+j-1;
    if(j<p)
        fitlist.phi_0(j,j+1:p,1,:) = beta(start+1:start+p-j,:);
        fitlist.phi_0(j+1:p, j,1,:) = beta(start+1:start+p-j,:);
    end
    start = start+p-j;
    for k = 1:q
        if(j>1)
            fitlist.phi_j(j,1:j-1,k,3,:) = beta(start+1:start+j-1,:);
            fitlist.phi_j(1:j-1,j,k,3,:) = beta(start+1:start+j-1,:);
        end
        start = start+j-1;
        if(j<p)
            fitlist.phi_j(j,j+1:p,k,2,:) = beta(start+1:start+p-j,:);
            fitlist.phi_j(j+1:p, j,k,2,:) = beta(start+1:start+p-j,:);
        end
        start = start+p-j;
    end
end
end

    