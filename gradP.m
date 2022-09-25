% Matlab function m-file gradP.m
%
% INPUT: 
%
% n-dimensional vector x
%
% OUTPUT:
%
% n-dimensional gradient vector of P evaluated at s

function grad = gradP(pk,s,X)
S=repmat(s(1:2),size(X,1),1);
D=S-X;
grad(1) = 2*sum(D(:,1)) + 2*pk*sum(D(:,1).*subplus((sum(D.^2,2)-s(3))));
grad(2) = 2*sum(D(:,2)) + 2*pk*sum(D(:,2).*subplus((sum(D.^2,2)-s(3))));
grad(3) =1-pk*sum(subplus(sum(D.^2,2)-s(3)));

end