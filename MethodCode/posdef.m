
function mnew = posdef(m,force)
% Function to add a multiple of identical matrix to make a given symmetric
% matrix positive definite. The multiplier is 2*min(eig(m)).
% Input: 
% m: a symmetric matrix but not necessarily positive definite
% Output: 
% mnew: a positive matrix with off diagonal elements same as m and diagonal
% elements differing a common factor.

p = size(m,1);
if m ~= m'
    error('input matrix is not symmetric');
elseif sum(sum(abs(m)))==0
    mnew = m;
elseif force == 0
   a = min(eig(m));
   if a>=0
       mnew = m;
   elseif a<0;
       mnew = m-5/4*a*eye(p); % if the min eigvalue is negative, make it positive by adding -2tims to all diagonal elements
   end
elseif force == 1
   a = min(eig(m));
   if a>0
       mnew = m;
   elseif a<0;
       mnew = m-5/4*a*eye(p); % if the min eigvalue is negative, make it positive by adding -2tims to all diagonal elements
   else
       mnew = m+abs(mean(mean(m)))/5*eye(p);
   end
   
%   k = 0;
%    while cond(mnew)>=10^4 % if the condition number is too large, add more to diagonal
%        mnew = mnew-2^k*a*eye(p);
%        k = k +1;
%    end
end
end

       
       
       
       
   
