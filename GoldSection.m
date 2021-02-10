% Optimization Theory
% Class Assignment  - BFGS
% 27/12/2020



function Return_Degeri =GoldSection(functname,tol,x,s,lowbound,intvl,ntrials)
format compact;

% find upper bound
upval = UpperBound(functname,x,s,lowbound,intvl,ntrials);
au=upval(1);	fau = upval(2);

% if upper bound returns value close to lowbound 
% return to the calling procedure to reverse the direction of the 
% search vector and try again

if (au <= 1.0e-06)
   aL = lowbound;  xL = x + aL*s;  
   faL =feval(functname,xL);
   Return_Degeri =[aL faL x];
   return
end


if (tol == 0) tol = 0.0001;  %default
end

eps1 = tol/(au - lowbound);
tau = 0.38197;
nmax = round(-2.078*log(eps1)); % no. of iterations




aL = lowbound;              xL = x + aL*s;  faL =feval(functname,xL);	
a1 = (1-tau)*aL + tau*au;   x1 = x + a1*s; fa1 = feval(functname,x1);
a2 = tau*aL + (1 - tau)*au; x2 = x + a2*s; fa2 = feval(functname,x2);

% storing all the four values for printing 
% remember to suppress printing after debugging
%fprintf('start  \n')
%fprintf('alphal(low)   alpha(1)   alpha(2)  alpha{up) \n')
avec = [aL a1 a2 au;faL fa1 fa2 fau];
%disp([avec])
for i = 1:nmax

	if fa1 >= fa2
   	aL = a1;	faL = fa1;
   	a1 = a2;	fa1 = fa2;
      a2 = tau*aL + (1 - tau)*au;	x2 = x + a2*s; 
          fa2 = feval(functname,x2);

   	au = au;	fau = fau;  % not necessary -just for clarity
      
      
      %fprintf('\niteration '),disp(i)
      %fprintf('alphal(low)   alpha(1)   alpha(2)  alpha{up) \n')
      avec = [aL a1 a2 au;faL fa1 fa2 fau];
		%disp([avec])

	else
  	   au = a2;	fau = fa2;
   	a2 = a1;	fa2 = fa1;
      a1 = (1-tau)*aL + tau*au;	x1 = x + a1*s; 
          fa1 = feval(functname,x1);
   	aL = aL;	faL = faL;  % not necessary
      
      %fprintf('\niteration '),disp(i)
      %fprintf('alphal(low)   alpha(1)   alpha(2)  alpha{up) \n')
      avec = [aL a1 a2 au;faL fa1 fa2 fau];
		%disp([avec])

	end
end

fprintf('alpha: ')
%avec = [aL a1 a2 au;faL fa1 fa2 fau];
avec=[a1];
disp([avec])

% returns the value at the last iteration
Return_Degeri =[a1 fa1 x1];







