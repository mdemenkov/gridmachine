function [Sol,Stat]=GridMachine(Func,FuncData,opt)

% [Sol,Stat]=GridMachine(Func,FuncData,Options)
%
% Return all solutions for a system of nonlinear equations inside given box
% Recommended for smooth equations only due to use of Newton's method
% Algorithm (c) by Mikhail Goman, Friday 13 of June 2003
% 
% Func - function handle or name, calling convention is Func(x,k,FuncData)
% where FuncData - any data, k-component number and x - vector of variables on which 
% the function depends
%
% The following parameters must be provided:
%
% Options.Limits(:,1)- minimal limits for x
% Options.Limits(:,2)- maximal limits for x
% Options.Density - density of a rectangular grid (the same for each
%                   coordinate axis)
% Options.DeriSteps - vector of derivative steps 
%
% The following parameters can be omitted:
%
% Options.SamePoints -tolerance of the same points 
% Options.NewtonZero - tolerance for zero in Newton's method
% Options.JacStep - re-compute Jacobian only after JacStep steps
% Options.MaxTime - number of maximum allowable iterations of Newton's
% Options.Progress==1 - display progress bar
% Options.ProgTitle - progress bar title
%
% Sol - each column of this output matrix represents one solution
%
% Stat.AllCubes - number of all generated cubes
% Stat.ProcessedCubes - number of all 'suspected' cubes
% Stat.RedundantSolutions - number of solutions which have been removed as
%                           redundant
% Stat.NewtonFailed - number of all cases when Newton algorithm was unable
%                     to find a solution

Limits=opt.Limits;
Density=opt.Density;

n=size(Limits,1);

Sol=[];
if isfield(opt,'SamePoints'), SamePoints=opt.SamePoints; else, SamePoints=1e-2; end
if isfield(opt,'Progress'), Progress=opt.Progress; else, Progress=1; end
if isfield(opt,'ProgTitle'), ProgTitle=opt.ProgTitle; else, ProgTitle='Computing...'; end

% Prepare options for Newton subroutine
DeriSteps=opt.DeriSteps;
if isfield(opt,'JacStep'), JacStep=opt.JacStep; else, JacStep=1; end
if isfield(opt,'MaxTime'), MaxTime=opt.MaxTime; else, MaxTime=50; end
if isfield(opt,'NewtonZero'), NewtonZero=opt.NewtonZero; else, NewtonZero=1e-5; end

% Obtain precomputed vertices of the n-cube
Vertices=GetCubeVertices(n);

NVertices=size(Vertices,1);
NCubes=Density^n;

cube=ones(n,1);dx=zeros(n,1);xgrid=zeros(n,Density+1);

for i=1:n
    dx(i)=(Limits(i,2)-Limits(i,1))/Density;
    xgrid(i,:)=Limits(i,1):dx(i):Limits(i,2);
end
% Progress indicator
if Progress
   hIBar = waitbar(0,ProgTitle);drawnow;
   binc=0.05;bcount=binc;
end

cubevert=zeros(NVertices,n);
x=zeros(n,1);

WarningState=warning;
warning off;

PCubes=0;RSols=0;NFailed=0;

for count=1:NCubes
    % Translate precomputed abstract cube vertices into real ones
    for i=1:NVertices
        for j=1:n
            if Vertices(i,j)
               cubevert(i,j)=xgrid(j,cube(j)+1);
            else
               cubevert(i,j)=xgrid(j,cube(j));
            end
        end
    end
    Outcast=0;
    for i=1:n
        range=feval(Func,cubevert',i,FuncData)';
        if sign(min(range))==sign(max(range))
           Outcast=1;break;
        end
    end
    % Outcasting cubes  that cannot contain intersection
    % of all nonlinear manifolds
    if ~Outcast
          PCubes=PCubes+1;
          cubecenter=(sum(cubevert)./NVertices)';
          [x,Flag]=Newton(Func,FuncData,cubecenter,DeriSteps,JacStep,MaxTime,NewtonZero);
          if Flag
                % Check if the solution is inside the given cube
                if all(Limits(:,1)<=x) && all(x<=Limits(:,2))
                   EFlag=0;
                   % Check if this is new solution 
                   nsol=size(Sol,2);
                   for j=1:nsol
                       if max(abs(Sol(:,j)-x))<=SamePoints
                          EFlag=1;RSols=RSols+1;break;
                       end
                   end
                   if ~EFlag, Sol=[Sol,x]; end
                end
          else
              NFailed=NFailed+1;
          end
    end % Outcasting cubes
    % Moving to the next cube
    for i=1:n
        if cube(i)<Density
           cube(i)=cube(i)+1; 
           break;
        else
           cube(i)=1;
        end
    end
    % Display progress
    if Progress
       if (count/NCubes)>=bcount
          waitbar(count/NCubes);bcount=bcount+binc;
       end
    end
end
if Progress, close(hIBar); end
warning(WarningState);

Stat.AllCubes=NCubes;
Stat.ProcessedCubes=PCubes;
Stat.RedundantSolutions=RSols;
Stat.NewtonFailed=NFailed;

function X=GetCubeVertices(n)

% Generate vertices of the n-cube
vert=zeros(1,n);bdigit=0;X=zeros(2^n,n);
for i=1:2^n
    for j=1:n
        if bitget(bdigit,j), vert(j)=1; else, vert(j)=0; end
    end
    X(i,:)=vert;
    bdigit=bdigit+1;
end

function [Sol,Flags]=Newton(Func,FuncData,x0,DeriSteps,JacStep,MaxTime,Tol)

% [x,flag]=Newton(Func,FuncData,x0,DeriSteps,JacStep,MaxTime,Tol)
% Classical Newton method
% 
% DeriSteps - derivative steps (different for each component of x)
% JacStep - re-compute Jacobian only after jacstep steps
% MaxTime - number of maximum allowable iterations
% Tol - x is the solution if all(abs(Func(x,FuncData))<=Tol)
%
% flag==1 if the solution x has been found

n=size(x0,1);m=size(x0,2);
Jac=zeros(n,n);Sol=zeros(n,m);Flags=zeros(1,m);
fx=zeros(n,1);

lastwarn('');
          
for i=1:m
    
    x=x0(:,i);
    for k=1:n
        fx(k)=feval(Func,x,k,FuncData);
    end
    iter=0;njac=JacStep;

    flag=0;JacIter=0;
    while any(abs(fx)>Tol) && iter<MaxTime
          if njac==JacStep
             njac=1;
             % Jacobian calculation
             for j=1:n
                 dxc=zeros(n,1);dxc(j)=DeriSteps(j);
                 for k=1:n
                     Jac(k,j)=(feval(Func,x+dxc,k,FuncData)-feval(Func,x-dxc,k,FuncData))./(2*DeriSteps(j));
                 end
             end
             JacIter=JacIter+1;
          else
             njac=njac+1;
          end
          % Obtain next point
          x_new=Jac\(Jac*x-fx);
          % Detect bad rank via 'warning'
          if isempty(lastwarn) % i.e. if rank(Jac)==n
             x=x_new;
          else
             break;
          end
          for k=1:n
              fx(k)=feval(Func,x,k,FuncData);
          end
          iter=iter+1;
    end
    flag=all(abs(fx)<=Tol);
    Sol(:,i)=x;Flags(i)=flag;
end
