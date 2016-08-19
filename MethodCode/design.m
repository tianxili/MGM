
function x = design(z, y, type, j)
% Function to obtain the design matrix for each separate regression
% Input: 
% z: data discrete
% y: data continuous
% type: 'z' or 'y'
% j: index of the regressed variable
% Output:
% x: the design matrix for the corresponding regression

n = size(z,1);
p = size(y,2);
q = size(z,2);
p1 = (p-1)*p/2;
if type == 'z'
    x = zeros(n,q-1+p+p1);
    tmp = z;
    tmp(:,j) = []; %% get z(-j)
    x(:,1:q-1) = tmp; %% the first q-1 columns are z(-j)
    x(:,q:q+p-1) = y;  %% next p columns are y
    ind1 = repmat(1:p,p,1);
    ind1 = ind1(tril(true(size(ind1)),-1)); %% if including y^2 terms, no -1
    ind2 = repmat(1:p,p,1)';
    ind2 = ind2(tril(true(size(ind2)),-1));
    yy = y(:,ind1).*y(:,ind2);
    x(:,q+p:q+p+p1-1) = yy;
elseif type == 'y'
    x = zeros(n,p+p*q-1);
    x(:,1:q) = z;
    tmp = y;
    tmp(:,j) = [];
    x(:,q+1:q+p-1) = tmp;
    x(:,q+p: end) = reshape(repmat(z, p-1,1),n,(p-1)*q).*repmat(tmp, 1,q);
end
clear tmp tmp1 ind1 ind2 yy
end
    
    
    
    
    
    
    
    
    
    
   