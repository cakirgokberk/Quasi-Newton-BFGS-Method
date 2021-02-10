
function RTV = myfunc(x)
% Optimization Theory
% Class Assignment  - BFGS
% 27/12/2020
% Objective Function for BFGS





% Baslangıc Noktasini Main.m 'de [-2 2] yapın. 
%RTV = 100*(x(2)-x(1).^2).^2 + (1-x(1)).^2           ;

% Baslangıc Noktasini Main.m 'de [1 2] yapın.
RTV = x(1).^2+0.5*x(2).^2+3                         ; % Quadratic

end
