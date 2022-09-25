
%
% INPUT: 
%
%   P              - the multivariable function to minimise (a separate
%                       user-defined Matlab function m-file)
%
%  gradP        - function which returns the gradient vector of f evaluated
%                       at x (also a separate user-defined Matlab function
%                       m-file)
%
% s0                - the starting iterate
%
% tolerance1   - tolerance for stopping criterion of steepest descent
%                           method
%
% tolerance2  - tolerance for stopping criterion of line minimisation (eg: in golden section search)
%
% T                     - parameter used by the "improved algorithm for
%                           finding an upper bound for the minimum" along
%                           each given descent direction
%
% OUTPUT:
%
% xminEstimate  - estimate of the minimum
%
% fminEstimate  - the value of f at xminEstimate

function [sminEstimate, PminEstimate,k] = BFGS(P, gradP, pk,s0,X,H0,tolerance1, tolerance2, T)
tic %starts timer

k = 0;  % initialize iteration counter
iteration_number=0; %initialise count

sk = s0;
sk_old=s0;
H_old=H0;
while ( norm(feval(gradP,pk, sk,X)) >= tolerance1 )
   iteration_number = iteration_number + 1;   
    % calculate steepest descent direction vector at current iterate xk
    %du = feval(gradf, xk);
    %dk = -du/norm(du);
    
    H_old = H_old / max(max(H_old)); % Correction if det H_old gets too lagre or small
    
    dk = transpose(-H_old*transpose(feval(gradP,pk, sk,X))); %gives dk as a row vector
    
    % minimise f with respect to t in the direction dk, which involves
    % two steps:
    
            % (1) find upper and lower bound,    [a,b]   ,for the stepsize t using the "improved
            % procedure" presented in the lecture notes
    
            [a, b] = multiVariableHalfOpen(P, pk, sk, X, dk, T);
        
               
             % (2) use golden section algorithm (suitably modified for
             % functions of more than one variable) to estimate the 
             % stepsize t in [a,b] which minimises f in the direction dk
             % starting at xk
                 
             [tmin, Pmin] = multiVariableGoldenSectionSearch(P, a, b, tolerance2, pk, sk, X, dk);
             
             % note: we do not actually need fmin, but we do need
             % tmin
             
             % update the steepest descent iteration counter and the
             % current iterate
             
             k = k + 1;
             
             sk = sk + tmin*dk;
             sk_new = sk_old +tmin*dk;
             
             csk=(sk_new -sk_old)'; %column vector
            
             gk= (feval(gradP,pk,sk_new,X)-feval(gradP,pk,sk_old,X))'; %column vector
             
             rk=(H_old*gk)/(csk'*gk);
             
             H_new = H_old + (1 + rk'*gk)/(csk'*gk)*(csk*csk') - (csk*rk') - (rk*csk');
             
             sk_old=sk_new;
             
             H_old=H_new;
             
end  
   
    % assign output values
    
   
    
    sminEstimate = sk;
    PminEstimate = feval(P,pk,sminEstimate,X);
