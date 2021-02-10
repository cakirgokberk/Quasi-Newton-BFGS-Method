
function RTV = myfunc(x)
% Optimization Theory
% Class Assignment  - BFGS
% 27/12/2020
% Objective Function for BFGS




%RTV = 100*(x(2)-x(1).^2).^2 + (1-x(1)).^2           ;

% Starting Point: [1,2] --> modifiable
RTV = x(1).^2+0.5*x(2).^2+3                         ; % Quadratic

end
