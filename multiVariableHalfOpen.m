% Matlab function m-file:  multiVariableHalfOpen.m
%
% INPUT: 
%
%   f           multivariable function to minimise
%
%   x           starting point
%
%   d           direction vector
%
%   T           upper bound increment parameter
%
% OUTPUT:
%
% lower and upper bound [a,b] on location of minimum of f in the direction d from
% the point x.

function [a,b] = multiVariableHalfOpen(P, pk, s, X, d, T)

 k = 1; 
 
      p = s;
      q = s + T*d;
    
      fp = feval(P,pk, p, X);
      fq = feval(P,pk, q, X);
    
            while ( fp > fq )
    
                  k = k +1;
            
                   p = q;
                   fp = fq;
            
                  q = p + (2^(k-1))*T*d;
                   fq = feval(P, pk, q, X);
                 
            end
            
            if ( k == 1)
                
                a = 0;
                b = T;
                
            elseif (k == 2)
                
                a = 0;
                b = 3*T;
                
            else    
          
             u = [0:k-1];
             v = [0:k-3];
             
             a = T*sum(2.^v);           
             b = T*sum(2.^u);
             
             end
       