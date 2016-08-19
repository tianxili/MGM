function [prob cparlist] = zprob(parlist)
% Function to generate the probability distribution of z as well as the canonical parameters.
% Input: 
% parlist: the parameter list
% Output:
% prob: the probability distribution of z
% cparlist: canonical parameter structure (g(2^q,1), h(p, 2^q), K(p,p,2^q));
[p q] = size(parlist.ita_j);
prob = zeros(2^q,1);
cparlist = struct('g',zeros(2^q,1),'h',zeros(p,2^q),'K', zeros(p,p,2^q));
z = de2bi(0:2^q-1);
for j = 1:2^q
    ztmp = z(j,:);
    zz = ztmp'*ztmp;
    zztmp = zz(tril((true(size(zz))),-1))';
    %g = ztmp*parlist.lambda_j + zztmp*parlist.lambda_jk;
    cparlist.g(j) = ztmp*parlist.lambda_j + zztmp*parlist.lambda_jk;
    %h = parlist.ita_j*ztmp'+parlist.ita_jk * zztmp';
    cparlist.h(:,j) = parlist.ita_0 + parlist.ita_j*ztmp';
    z3 = reshape(kron(ztmp, ones(p,p)),[p p q]);
    %zz3 = reshape(kron(zztmp, ones(p,p)), [p p q*(q-1)/2]);
    %K = parlist.phi_0+ sum(z3.*parlist.phi_j,3)+sum(zz3.*parlist.phi_jk,3);
    cparlist.K(:,:,j) = parlist.phi_0+ sum(z3.*parlist.phi_j,3);
    prob(j) = det(cparlist.K(:,:,j))^(-0.5)*exp(cparlist.g(j)+cparlist.h(:,j)'*inv(cparlist.K(:,:,j))*cparlist.h(:,j)/2);
end
prob = prob/sum(prob);
end

    
    
    
    
    
    

