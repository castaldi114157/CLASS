function x = bisection(f, N, a, b)
%    Synopsis:
%    Rudimentary implementation of the Bisection Method. A function handle
%    is used to evaluate the objective function f(x) whose root(s) you want
%    to find. 
%
%    inputs:
%    f            f(x)
%    N            max number of iterations
%    a            lower bound
%    b            upper bound
%
%    output:      
%    x            the value of the root if the method is successful.


% initialize the iteration counter
k = 0;
% set a default tolerance
tol = 1.0e-6;
% set a dummy value for the solution x  
x = 0; 


% Sane Checks:
% ============
if (a > b)
    % check if the user has inverted the lower and upper bounds
    % and correct it if that's the case
    tmp=a;
    a = b;
    b = tmp;
end

% verify that there's a indeed a root in the given interval    
if ( f(a) * f(b) > 0 )
    fprintf('>> No roots found in the given interval [a, b].\n');
    fprintf('>> Change [a, b] and try again.');
    fprintf('>> Bisection Method Terminated.\n');
    return    
end

% ----------------------------------------------------------------------- %
% verify if the user has already provided a satisfactory solution that 
% meets the tolerance specifications. 
if (abs( f(a) ) < tol) 
    % x = a is a root to the speicified tolerance
    x = a; 
    fprintf('>> Solution found in %d iterations.\n', k);
    fprintf('>> x = %g \t\t f(x) = %g\n', x, f(a));
    return
end

if (abs( f(b) ) < tol) 
    % x = b is a root to the specified tolerance
    x = b; 
    fprintf('>> Solution found in %d iterations.\n', k);
    fprintf('>> x = %g \t\t f(x) = %g\n', x, f(b));
    return
end
% ----------------------------------------------------------------------- %

% status; used to check if the Method failed to find a solution with the
% gfiven number of iterations (N).
status=1; 
for k = 1:N
    % Bisection Method
    
    % estimate the root with the middle value and evaluate f(x)
    x = (a * f(b) - b * f(a))/(f(b) - f(a));
    % interval length
    L = b - a;
    % display iterations
    fprintf('i = %2d \t\t a = %+.4e,  f(a) = %+.4e \t\t b = %+.8f,  f(b) = %+.4e \t\t x = %+.4e,  f(x) = %+.4e \t\t  L = %+.4e \n',...
        k, a, f(a), b, f(b), x, f(x), L);
    
    % check if the root estimate meets the convergence criteria 
    if ( abs(f(x)) <= tol || L <= tol )
        % notify user that a solution has been found
        fprintf('>> Solution found in %d iterations\n', k);
        fprintf('>> x = %g \t f(x) = %g\n', x, f(x));
        status=0;
        break;
    end
    
    % assume that x is the new upper bound
    if ( f(a) * f(x) < 0 )
        % set the new upper bound
        b = x; 
    else
        % otherwise x is the lower bound
        a = x;
    end
end

if (status == 1)
    % notify user if a solution could not be found within the max number of
    % iterations
    fprintf(['>> Bisection Method failed to find ',...
       'a solution in %d iterations'], N); 
end

return