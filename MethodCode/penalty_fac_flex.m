%% function to get the weight factor for the covariates in a regression. 
% see paper for the definitions of weights. 

function w = penalty_fac_flex(type,p,q,kappa)
p1 = (p-1)*p/2;
if strcmp(type, 'z')
    w = zeros(p+q+p1-1,1);
    w(1:(q-1)) = 1*kappa;
    w(q:q+p-1) = 1;
    w(q+p:end) = 2;
elseif strcmp(type, 'y')
    w = zeros(p-1+p*q,1);
    w(1:q) = 1;
    w(q+1:q+p-1) = 1;
    w(q+p:end) = 2;
end