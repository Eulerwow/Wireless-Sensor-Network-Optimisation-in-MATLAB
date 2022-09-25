function val = f(s,X)
S=repmat(s(1:2),size(X,1),1);
D=S-X;
val= sum(sum(D.^2)) + max(sum(D.^2,2));
end 