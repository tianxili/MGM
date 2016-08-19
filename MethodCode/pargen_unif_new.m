
function parlist = pargen_unif_new(adjmat, p, q, a, b, c)
%% Get the sub adjascent matrix structure for zz, zy, yy;
if (p+q) ~= size(adjmat,1)
    error('dimension dismatch');
end
adj = struct('zz',{},'zy',{},'yy',{});
adj(1).zz = adjmat(1:q,1:q);
adj(1).zy = adjmat(1:q, (q+1):(q+p));
adj(1).yy = adjmat((q+1):(q+p), (q+1):(q+p));

%% Generate the intial random parameters
parlist = struct('lambda_j',{},'lambda_jk',{},'ita_0',{},'ita_j',{}, 'phi_0', {}, 'phi_j', {});
parlist(1).lambda_j = c*sign(unifrnd(-1,1,q,1)).*unifrnd(0.9,1.1,q,1);
parlist(1).lambda_jk = c*sign(unifrnd(-1,1,q*(q-1)/2,1)).*unifrnd(0.9,1.1,q*(q-1)/2,1);
parlist(1).ita_0 = a*sign(unifrnd(-1,1,p,1)).*unifrnd(0.9,1.1,p,1);
parlist(1).ita_j = a*sign(unifrnd(-1,1,p,q)).*unifrnd(0.9,1.1,p,q);
%tmp = rand(p);
%tmp = triu(unifrnd(a/2,a,p,p).*sign(unifrnd(-1,1,p,p)));
tmp = triu(b*unifrnd(0.9,1.1,p,p),1);
parlist(1).phi_0 = tmp + tmp';
parlist(1).phi_j = zeros(p,p,q);
for j = 1:q
    %tmp = rand(p);
    %tmp = triu(unifrnd(a/2,a,p,p).*sign(unifrnd(-1,1,p,p)));
    tmp = triu(b*unifrnd(0.9,1.1,p,p),1);
    parlist(1).phi_j(:,:,j) = tmp + tmp';
end

%% Reset some of the parameters to 0 according to the adjascent matrix;
%%% edges of Zj and Zk
for j = 1:(q-1)
    for k = (j+1):q
        if(adj.zz(j,k)==0)
            ind = transind(q,j,k);
            parlist.lambda_jk(ind) = 0;
        end
    end
end
%%% edges of Zj and Yk
for j = 1:q
    for k = 1:p
        if adj.zy(j,k)==0
        parlist.ita_j(k,j) = 0;
        parlist.phi_j(:, k, j) = 0;
        parlist.phi_j(k, :, j) = 0;
        end
    end
end
%%% edges of Yj and Yk
for j = 1:(p-1)
    for k = (j+1):p
        if adj.yy(j,k)==0
            parlist.phi_0(j,k) = 0;
            parlist.phi_0(k,j) = 0;
            parlist.phi_j(j,k,:) = 0;
            parlist.phi_j(k,j,:) = 0;
        end
    end
end
% Readjust the Phi matrices to make them positive definite;
parlist.phi_0 = posdef(parlist.phi_0,1);
for j = 1:q
    parlist.phi_j(:,:,j) = posdef(parlist.phi_j(:,:,j),0);
    parlist.phi_0 = parlist.phi_0 + diag(diag(parlist.phi_j(:,:,j)));
    parlist.phi_j(:,:,j) = parlist.phi_j(:,:,j) - diag(diag(parlist.phi_j(:,:,j)));
end

% parlist.lambda_j = abs(parlist.lambda_j);
% parlist.lambda_jk = abs(parlist.lambda_jk);
% parlist.ita_0 = abs(parlist.ita_0);
% parlist.ita_j = abs(parlist.ita_j);

end



            


            






                
