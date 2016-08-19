function [z y prob cparlist] = datagen(parlist, n)
% Function to generate data given the parameters
[p q] = size(parlist.ita_j);
[prob cparlist] = zprob(parlist);
a = 0;
b = 0;
while (a==0 || b==1)
ztmp = randp(prob,n,1);            
z = de2bi(ztmp-1, q);
sprintf('binary z generation finished')
a = min(mean(z,1));
b = max(mean(z,1));
end
y = zeros(n,p);
for i = 1:n
    sigma = inv(cparlist.K(:,:,ztmp(i)));
    mu = cparlist.K(:,:,ztmp(i))\ cparlist.h(:,ztmp(i));
    y(i,:) = mvnrnd(mu, sigma);
end
sprintf('data generation finished')
end

