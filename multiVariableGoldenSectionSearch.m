% Matlab function m-file : multiVariableGoldenSectionSearch.m
%
% performs golden section search for finding minimum of f along the
% direction d, starting at x, where the minimum has upper and lower bound
% [a,b]


function [minEstimate, PminEstimate] = multiVariableGoldenSectionSearch(P,a,b,tolerance,  pk, s, X, d)

% Check that the correct number of input arguments have been passed to the
% function (this is not essential, but it helps to detect user error)

if (nargin ~=8)     % in Matlab the operator ~= means "not equal to"
    error('six input arguments are required')
end

% Check that input parameters have appropriate values

    if(b <= a)
        error('b must be strictly greater than a')    
    end
    
    if(tolerance <= 0)
        error('tolerance must be strictly positive')
    end

    
% begin the Golden Search algorithm

gamma = (sqrt(5) - 1)/2;  

k = 1;  % iteration counter

p = b - gamma*(b-a);
q = a + gamma*(b-a);

fp = feval(P, pk, s + p*d, X);
fq = feval(P, pk, s + q*d, X);

while ( (b-a) >= 2*tolerance )
    
    k = k + 1;
    
    if (fp <= fq) 
        
        b = q;
        q = p;
        fq = fp;
        
        p = b - gamma*(b-a);
        fp = feval(P, pk, s + p*d, X);
        
    else   
    
        a = p;
        p = q;
        fp = fq;
        
        q = a + gamma*(b-a);
        fq = feval(P, pk, s + q*d, X);
        
    end

end

% assign the output values of this function

minEstimate = (a+b)/2;      % the midpoint of the final interval
PminEstimate = feval(P, pk, s + minEstimate*d, X);  


    
    
    
    
    
    
    
    
    
    
    