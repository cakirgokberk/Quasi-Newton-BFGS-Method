function ReturnV = UpperBound(functname,x,s,a0,da,ns)
% Optimization Theory
% Class Assignment  - BFGS
% 27/12/2020
% the current position vector 			x
% the current search direction vector	s
% the initial step						a0
% the incremental step					da
% the number of bracketting steps		ns

format compact
%	ntrials are used to bisect/double values of da
if (ns ~= 0) ntrials = ns;
else ntrials = 10;   % default
end

if (da ~= 0) das = da;
else das = 1;  %default
end
% finds a value of function greater than or equal
% to the previous lower value

for i = 1:ntrials
   j = 0;	dela = j*das;	a00 = a0 + dela;  
   dx0 = a00*s;	x0 = x + dx0;  f0 = feval(functname,x0);
   j = j+1;	dela = j*das;	a01 = a0 + dela;
   dx1 = a01*s;	x1 = x + dx1;	f1 = feval(functname,x1);
   f1s = f1;
   if f1 < f0 
         for j = 2:ntrials
         	a01 = a0 + j*das;		dx1 = a01*s;	
            x1 = x + dx1;		f1 = feval(functname,x1);
            f1s = min(f1s,f1);
            if f1 > f1s 
      			ReturnV = [a01 f1 x1];
               return;
            end
         end
         %fprintf('\nCannot increase function in ntrials')
      	ReturnV = [a01 f1 x1];
         return;
			         
   elseif	f1 >= f0
      das = 0.5*das;
   end
end
ReturnV =[a0 f0 x0];
   
   