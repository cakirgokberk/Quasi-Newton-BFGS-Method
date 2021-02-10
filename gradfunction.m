function Return = gradfunction(functname,x)
% Optimization Theory
% 18/12/2020

hstep = 0.001;
n = length(x);
f = feval(functname,x);

for i = 1:n
   xs = x;
   xs(i) = xs(i) + hstep;
   gradx(i)= (feval(functname,xs) -f)/hstep;
end
Return = gradx;