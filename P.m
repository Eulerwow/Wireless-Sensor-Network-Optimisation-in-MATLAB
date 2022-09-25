%Input s: 3x1 vector ,pk: natural number, and X nx2 matrix 
function val = P(pk,s,X)
S=repmat(s(1:2),size(X,1),1);
D=S-X;
val= sum(sum(D.^2)) + s(3)+(pk/2)*sum( ...
    ( subplus(    sum(D.^2,2)-s(3)    ) ).^2 ...
    );
end  