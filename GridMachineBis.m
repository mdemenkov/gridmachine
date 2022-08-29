function [Sol,Stat]=GridMachineBis(Func,FuncData,opt)

% [Sol,Stat]=GridMachineBis(Func,FuncData,Options)
%
% Return all solutions for a system of nonlinear equations 
% Bisection method, well suited for nonconvex equations
% Algorithm (c) by Mikhail Goman, Monday 16 of June 2003
% 
% Func - function handle or name, calling convention is Func(x,k,FuncData)
% where FuncData - any data, k-component number and x - vector of variables on which 
% the function depends
%
% The following parameters must be provided:
%
% Options.Limits(:,1)- minimal limits for x
% Options.Limits(:,2)- maximal limits for x
% Options.StackSize - maximal number of generated cubes
% Options.InitialGrid - initial density of a rectangular grid (the same for each
%                       coordinate axis)
%
% The following parameters can be omitted:
%
% Options.ZeroCube - tolerance for zero
% Options.Progress==1 - display progress bar
%
% Sol - each column of this output matrix represents one solution
%
% Stat.MaxCubesInStack - peak number of generated cubes
% Stat.NumberOfSolutions - self-explained
% Stat.FinalCubes - number of cubes in the final stack at the end
% Stat.StackOverflow =1 if stack overflow occurs
% Stat.RedundantSolutions - number of solutions which have been removed as
%                           redundant
% Stat.Iterations - number of all bisections applied
% Stat.AllCubes - number of all generated cubes

n=size(opt.Limits,1);

StackSize=opt.StackSize;

CubeLim=opt.Limits;
SubCube=zeros(n,2);

CubeStack=zeros(StackSize*n,2);
StackMirror=CubeStack;
FinalStack=CubeStack;

StackCount=0;
MirrorCount=0;
FinalCount=0;

MaxCubes=0;RSols=0;Iteration=0;PCubes=0;

Sol=[];
if isfield(opt,'ZeroCube'), ZeroCube=opt.ZeroCube; else, ZeroCube=1e-5; end
if isfield(opt,'Progress'), Progress=opt.Progress; else, Progress=1; end
if isfield(opt,'SamePoints'), SamePoints=opt.SamePoints; else, SamePoints=1e-2; end
if isfield(opt,'InitialGrid'), InitialGrid=opt.InitialGrid; else, InitialGrid=1; end

% Obtain precomputed vertices of the n-cube
Vertices=GetCubeVertices(n);

NVertices=size(Vertices,2);

WarningState=warning;
warning off;

Sides=CubeLim(:,2)-CubeLim(:,1);
MaxIter=n*ceil(log2(max(Sides)/InitialGrid/ZeroCube));

disp(['Num. of iterations <= ' num2str(MaxIter)]);

%try
    
% ------------------------     Generate initial partition -----------------------

      dx=zeros(n,1);xgrid=zeros(n,InitialGrid+1);
      for i=1:n
          dx(i)=(CubeLim(i,2)-CubeLim(i,1))/InitialGrid;
          xgrid(i,:)=CubeLim(i,1):dx(i):CubeLim(i,2);
      end
      
      subcube=ones(n,1);
    
      NCubes=InitialGrid^n;
      
      for count=1:NCubes % Divide initial cube
          
          % Obtain coordinate limits for the current cube
          for i=1:n
              SubCubeLim(i,1)=xgrid(i,subcube(i));
              SubCubeLim(i,2)=xgrid(i,subcube(i)+1);
          end
          
          % Put this cube into stack
          StackCount=StackCount+1;
          if StackCount>StackSize
             error('Stack overflow');
          end
          CubeStack((StackCount-1)*n+1:StackCount*n,:)=SubCubeLim;
          if MaxCubes<StackCount, MaxCubes=StackCount; end
          
          % Moving to the next subcube
          for i=1:n
              if subcube(i)<InitialGrid
                 subcube(i)=subcube(i)+1; 
                 break;
              else
                 subcube(i)=1;
              end
          end
          
      end % Divide cube

StackCount=NCubes;

% -------------------------- Main cycle --------------------------------------


while StackCount>0  % Main cycle
 
  PCubes=PCubes+StackCount;
  
  MirrorCount=0;Iteration=Iteration+1;

  if Progress
     ProgTitle=['Iteration=' num2str(Iteration)...
                ',main stack=' num2str(StackCount) ' cube(s)'...
                ', final stack=' num2str(FinalCount) ' cube(s)'...
                ', max. size=' num2str(MaxCubeSize(CubeStack,StackCount,n))]; 
     disp(ProgTitle);drawnow;
  end
      
  Volume=0;
  
  for StackInd=1:StackCount % Stack processing
      
      % Obtain coordinate limits for the current cube
      CubeLim=CubeStack((StackInd-1)*n+1:StackInd*n,:);
      
      % Split this cube along biggest side
      Sides=CubeLim(:,2)-CubeLim(:,1);
      ind=find(Sides==max(Sides));ind=ind(1);
      meanval=CubeLim(ind,1)+Sides(ind)/2;
      
      SubCube1=CubeLim;SubCube2=CubeLim;
      SubCube1(ind,2)=meanval;
      SubCube2(ind,1)=meanval;
      
      CubeVert1=zeros(n,NVertices);CubeVert2=CubeVert1;
      
      % Translate precomputed abstract cube vertices into real ones
      for i=1:NVertices
          for j=1:n
              if Vertices(j,i)
                 CubeVert1(j,i)=SubCube1(j,2);
                 CubeVert2(j,i)=SubCube2(j,2);
              else
                 CubeVert1(j,i)=SubCube1(j,1);
                 CubeVert2(j,i)=SubCube2(j,1);
              end
          end
      end
      
      for count=1:2 % Check subcubes
          
          if count<2
             SubCubeLim=SubCube1;
             CubeVert=CubeVert1;
          else
             SubCubeLim=SubCube2;
             CubeVert=CubeVert2;
          end
          
          % Outcasting cubes  that cannot contain intersection
          % of all nonlinear manifolds
          Outcast=0;
          for i=1:n
              range=feval(Func,CubeVert,i,FuncData);
              if sign(min(range))==sign(max(range))
                 Outcast=1;break;
              end
          end
          if ~Outcast
             Sides=SubCubeLim(:,2)-SubCubeLim(:,1);
             % Is it small enough ?
             if all(Sides<=ZeroCube)
                % Put this cube into final stack
                FinalCount=FinalCount+1;
                if FinalCount>StackSize
                   error('Final stack overflow');
                end
                FinalStack((FinalCount-1)*n+1:FinalCount*n,:)=SubCubeLim;
             else
                % Put this cube into main stack for futher processing
                MirrorCount=MirrorCount+1;
                if MirrorCount>StackSize
                   error('Main stack overflow');
                end
                StackMirror((MirrorCount-1)*n+1:MirrorCount*n,:)=SubCubeLim;
                if MaxCubes<MirrorCount, MaxCubes=MirrorCount; end
             end
          end 
         
      end % Check subcubes
            
  end  % Stack processing
  
  CubeStack(1:MirrorCount*n,:)=StackMirror(1:MirrorCount*n,:);
  StackCount=MirrorCount;
  
 figure;hold on;axis([0 10 0 10]);
 ShowStack(CubeStack,StackCount,n,1,2);
  
end 

% Cubes in the final stack should contain all solutions
if FinalCount>0
        NSol=0;RSols=0;
        for StackInd=1:FinalCount
            % Obtain coordinate limits for the current cube
            CubeLim=FinalStack((StackInd-1)*n+1:StackInd*n,:);
            Sides=CubeLim(:,2)-CubeLim(:,1);
            % Return its center as a solution
            x=CubeLim(:,1)+Sides./2;
            % Check if this is new solution
            EFlag=0;
            for i=1:NSol
                if max(abs(Sol(:,i)-x))<=SamePoints
                   EFlag=1;RSols=RSols+1;
                   break;
               end
            end
            if ~EFlag, Sol=[Sol,x];NSol=NSol+1; end
        end
end
 
%catch % Stack  overflow
    
 %   Stat.StackOverflow=1;
 %   warning(WarningState);
    
 %   return
    
 %end

warning(WarningState);

Stat.StackOverflow=0;
Stat.Iterations=Iteration;
Stat.NumberOfSolutions=NSol;
Stat.MaxCubesInStack=MaxCubes;
Stat.FinalCubes=FinalCount;
Stat.RedundantSolutions=RSols;
Stat.AllCubes=PCubes;

function X=GetCubeVertices(n)

% Generate vertices of the n-cube
vert=zeros(n,1);bdigit=0;X=zeros(n,2^n);
for i=1:2^n
    for j=1:n
        if bitget(bdigit,j), vert(j)=1; else, vert(j)=0; end
    end
    X(:,i)=vert;
    bdigit=bdigit+1;
end

function ShowStack(CubeStack,StackCount,n,i,j)

for StackInd=1:StackCount
            
    CubeLim=CubeStack((StackInd-1)*n+1:StackInd*n,:);
    plot(CubeLim(i,:),CubeLim(j,1)*[1 1],'b');
    plot(CubeLim(i,:),CubeLim(j,2)*[1 1],'b');
    plot(CubeLim(i,1)*[1 1],CubeLim(j,:),'b');
    plot(CubeLim(i,2)*[1 1],CubeLim(j,:),'b');
end

function d=CDist(Cube,x)

n=length(x);
int_dist=zeros(n,1);
for i=1:n
    if x(i)>Cube(i,2)
       int_dist(i)=x(i)-Cube(i,2);
    elseif x(i)<Cube(i,1)
       int_dist(i)=Cube(i,1)-x(i);
    else
       int_dist(i)=0;
    end
end
d=max(int_dist);

function d=CSize(Cube)

Sides=Cube(:,2)-Cube(:,1);
d=max(Sides);

function d=MaxCubeSize(CubeStack,StackCount,n)

f=zeros(StackCount,1);
for StackInd=1:StackCount
    CubeLim=CubeStack((StackInd-1)*n+1:StackInd*n,:);
    Sides=CubeLim(:,2)-CubeLim(:,1);
    f(StackInd)=max(Sides);
end
d=max(f);
