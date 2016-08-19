function l = transind(dim,j,k)

% Function to calculate the index of (j,k(i)) in the rearranged triangular
% array. k can be a vector. 
% Input: 
% dim: the total dimension of the original array
% j: first index
% k: second index, can be a vector
% Output:
% l: the index of (j,k) in the rearranged array, same length as k

l  = zeros(length(k),1);
if (j>dim) || (max(k)>dim)
    error('Incorrect Input')
else
    for i = 1:length(k)
        if j<k(i)
            l(i) = (j-1)*(2*dim-j)/2+k(i)-j;
        elseif j>k(i)
            l(i) = (k(i)-1)*(2*dim-k(i))/2+j-k(i);
        end 
    end
end

end
