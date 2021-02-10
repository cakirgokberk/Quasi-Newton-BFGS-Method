function Return = BFGS(functname,dvar0,niter,tol,lowbound,intvl,ntrials)
% Optimization Theory
% BFGS
% 27/12/2020


clc    
clf   
warning off  

e1 = 1.0e-04;           % convergence
e2 = 1.0e-08; 
e3 = 1.0e-04;  
nvar = length(dvar0);   % number of variables



%% Plotting

if (nvar == 2)
   
    delx1 = 6;
    delx2 = 5;

    x1 = (dvar0(1)-delx1):0.1:(dvar0(1)+delx1);
    x2 = (dvar0(2)-delx2):0.1:(dvar0(2)+delx2);
    x1len = length(x1);
    x2len = length(x2);
    
    for i = 1:x1len
        for j = 1:x2len
            x1x2 =[x1(i) x2(j)];
            fun(j,i) = feval(functname,x1x2);
        end
    end
   
    
    c1 = contour(x1,x2,fun,[3.1 3.25 3.5 4 6 10 15 20 25],'k');
      
    grid
    xlabel('x_1');
    ylabel('x_2');
    funname = strrep(functname,'_','-');
    title(strcat('BFGS:',funname));

  end


xs(1,:) = dvar0;
x = dvar0;
Lc = 'r';
fs(1) = feval(functname,x);             % Initial Value of Function
as(1)=0;
grad = (gradfunction(functname,x)); 

H = eye(nvar);                          % initial Q
convg(1)=grad*grad';




for i = 1:niter-1
    
    fprintf('iteration number:  '),disp(i)
    
    d = (-inv(H)*grad')';                                                   % Search Direction
    
    output = GoldSection(functname,tol,x,d,lowbound,intvl,ntrials);         % Finding Alpha
    as(i+1) = output(1);
    fs(i+1) = output(2);
    
    for k = 1:nvar
        xs(i+1,k)=output(2+k);
        x(k)=output(2+k);
    end
 
    
    grad= (gradfunction(functname,x)); 
    convg(i+1)=grad*grad';
    
    
    
    fprintf('objective function value:  '),disp(fs(i+1));
    
    
   
    %% Plotting (drawing lines)
    
    if (nvar == 2)
        line([xs(i,1) xs(i+1,1)],[xs(i,2) xs(i+1,2)],'LineWidth',2,'Color',Lc)
        itr = int2str(i);
        x1loc = 0.5*(xs(i,1)+xs(i+1,1));
        x2loc = 0.5*(xs(i,2)+xs(i+1,2));
        
       
        if strcmp(Lc,'r') 
            Lc = 'b';
        else
            Lc = 'r';
        end
        
        pause(0.5)  
       
    end
    
   %% METHOD 
    
    if(convg(i+1)<= e3)     % convergence criteria
        break; 
    
    
    end               
    
    
    
    delx = (x - xs(i,:))';                      % update the metric here
    gradold = gradfunction(functname,xs(i,:));
    Y = (grad -gradold)'; 
    B = (Y*Y')/(Y'*delx);
    C = gradold'*gradold/(gradold*d');
    H = H + B + C;
    
    
   
    
end


len=length(as);
designvar=xs(length(as),:);
fprintf('\n Iterations:  '),disp(len-1)
fprintf('\n     x(1)      x(2)      F(x)\n')
disp([xs fs']);
Return = [designvar fs(len)];































