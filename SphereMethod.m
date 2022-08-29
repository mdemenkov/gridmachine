function Branches=SphereMethod(Func,FuncData,Options)

% Branches=SphereMethod(Func,FuncData,Options)
%
% Func - function handle or name, calling convention is \dot x_k=Func(x,k,FuncData,cp)
% where FuncData - any data, k-component number, x - vector of variables on which 
% the function depends and cp - the continuation parameter
%
% The following parameters must be provided:
%
% Options.SolverFun - solver for nonlinear equations (function handle)
% Options.SolverOpt - parameters for the solver
%
% Options.StackSize - maximal number of simultaneously investigated branches
% Options.Box(:,1)- minimal limits for [x;cp]
% Options.Box(:,2)- maximal limits for [x;cp]
% Options.xcp0 - initial value for [x;cp]
% Options.MaxRad - maximum radius of the ball
% Options.MinRad - minimum radius of the ball
%
% The following parameters can be omitted:
%
% Options.SamePoints - tolerance of the same points 
% Options.ScaleCube - scale coefficient for surrounding cube
% Options.ScaleSphere - scale coefficient for reducing the sphere
% Options.Progress=1 - display information about computational process
%
% Output:
%
% Branches - cell array of branches as matrices where each column represent
%            a point

if isfield(Options,'SamePoints'), SamePoints=Options.SamePoints; else, SamePoints=1e-2; end
if isfield(Options,'ScaleCube'), ScaleCube=Options.ScaleCube; else, ScaleCube=1; end
if isfield(Options,'ScaleSphere'), ScaleSphere=Options.ScaleSphere; else, ScaleSphere=1/2; end
if isfield(Options,'Progress'), dflag=Options.Progress; else, dflag=1; end

Branches=[];
Box=Options.Box;

n=size(Box,1);

StackSize=Options.StackSize;
Stack=zeros(StackSize*n,2);
StackMirror=Stack;
StackCount=0;MirrorCount=0;BranchCount=0;

MaxRad=Options.MaxRad;
MinRad=Options.MinRad;

% Prepare parameters for ExtendedSystem
SystemData.Func=Func;
SystemData.FuncData=FuncData;

SolverFun=Options.SolverFun;
SolverOpt=Options.SolverOpt;

ex=Options.xcp0;

ddisp(dflag,'[x0;cp0]=');
ddisp(dflag,ex);

Sphere.Radius=MaxRad;
Sphere.Center=ex;
InitPoints=SphereStep(SolverFun,SolverOpt,SystemData,Sphere,[],ScaleCube,SamePoints);

if isempty(InitPoints)
   ddisp(dflag,'No initial points detected');
   return; 
else
   ddisp(dflag,'Initial points =');
   ddisp(dflag,InitPoints);
end

for i=1:size(InitPoints,2)
    StackCount=StackCount+1;
    Stack((StackCount-1)*n+1:StackCount*n,:)=[ex,InitPoints(:,i)];
end

% Prepare parameters for TraceOneBranch
TraceOpt.SolverFun=SolverFun;
TraceOpt.SolverOpt=SolverOpt;
TraceOpt.SystemData=SystemData;
TraceOpt.MaxRad=MaxRad;
TraceOpt.MinRad=MinRad;
TraceOpt.Box=Box;
TraceOpt.SamePoints=SamePoints;
TraceOpt.ScaleCube=ScaleCube;
TraceOpt.dflag=dflag;
TraceOpt.ScaleSphere=ScaleSphere;

% Main cycle
while StackCount>0
      ddisp(dflag,['Stack size=' num2str(StackCount)]);
      MirrorCount=0;
      for StackInd=1:StackCount
          StartPoints=Stack((StackInd-1)*n+1:StackInd*n,:);
          
          ddisp(dflag,'Starting points for a branch=');
          ddisp(dflag,StartPoints);
          ddisp(dflag,'Tracing the branch...');
                    
          [Sol,Stat]=TraceOneBranch(StartPoints,TraceOpt);
          
          if Stat.BranchingPoint
             ddisp(dflag,'Branching point detected');
          elseif Stat.EndOfBranch
             ddisp(dflag,'Branch disappeared');
          elseif Stat.OutOfBox
             ddisp(dflag,'Branch reached the limits');
          end
          
          if Stat.BranchingPoint
             ex=Sol.EndPoints(:,1);
             InitPoints=Sol.EndPoints(:,2:end);
             for i=1:size(InitPoints,2)
                 MirrorCount=MirrorCount+1;
                 if MirrorCount>StackSize
                    error('Stack overflow');
                 end
                 StackMirror((MirrorCount-1)*n+1:MirrorCount*n,:)=[ex,InitPoints(:,i)];
             end
          end
          if ~isempty(Sol.Branch)
             BranchCount=BranchCount+1;
             Branches{BranchCount}=Sol.Branch;
             ddisp(dflag,['Branch No ' num2str(BranchCount) ' has been collected']);
          end
      end
      Stack(1:MirrorCount*n,:)=StackMirror(1:MirrorCount*n,:);
      StackCount=MirrorCount;
end

function [Sol,Stat]=TraceOneBranch(StartPoints,TraceOpt)
%
% StartPoints - two starting points for the branch (first column is the
%               center of a sphere and second column is intersections
%               of the sphere with the branch)
%
% Sol.Branch - every column is a point of the branch
% Sol.EndPoints - starting points of new branches if branching point
%                 detected (first column is the center of a sphere and
%                 other points are intersections of the sphere with branches)
%
% Stat.BranchingPoint=1 - branching point detected
% Stat.EndOfBranch=1 - branch suddenly disappeared
% Stat.OutOfBox=1 - branch exceeded the limits

SolverFun=TraceOpt.SolverFun;
SolverOpt=TraceOpt.SolverOpt;
SystemData=TraceOpt.SystemData;
MaxRad=TraceOpt.MaxRad;
MinRad=TraceOpt.MinRad;
Box=TraceOpt.Box;
SamePoints=TraceOpt.SamePoints;
ScaleCube=TraceOpt.ScaleCube;
dflag=TraceOpt.dflag;
ScaleSphere=TraceOpt.ScaleSphere;

BranchingPoint=0;EndOfBranch=0;OutOfBox=0;
EndPoints=[];

nex=size(StartPoints,1);
Radius=MaxRad;
ex0=StartPoints(:,1);ex1=StartPoints(:,2);
Branch=ex1;

while ~(BranchingPoint | EndOfBranch | OutOfBox)

         dx=ex1-ex0;ex=ex1+nrm(dx)*Radius;
         Sphere.Center=ex;
         Sphere.Radius=Radius;
         Sol=SphereStep(SolverFun,SolverOpt,SystemData,Sphere,ex1,ScaleCube,SamePoints);
         if size(Sol,2)>1 % More than one new solution
             if Radius>MinRad
             % Contract the sphere to avoid capturing neighbouring branches
                Radius=Radius*ScaleSphere;
                if Radius<MinRad, Radius=MinRad; end
                ddisp(dflag,'Contracting sphere, radius=%0.5g',Radius);
             else
             % True branching point - expand the sphere to capture
             % emanating branches
                 np=size(Branch,2);
                 for i=np:-1:1
                     ex1=Branch(:,i);
                     dex=ex-ex1;
                     Radius=sqrt(dex'*dex);
                     if Radius>=MaxRad, break; end
                 end
                 ddisp(dflag,'Expanding sphere, radius=%0.5g',Radius);
                 Sphere.Center=ex;
                 Sphere.Radius=Radius;
                 Sol=SphereStep(SolverFun,SolverOpt,SystemData,Sphere,ex1,ScaleCube,SamePoints);
                 if size(Sol,2)>1
                    EndPoints=[ex,Sol];
                    BranchingPoint=1; break; 
                 else
                    error('Cannot locate branching point');
                 end   
             end
         elseif size(Sol,2)<1 % No solutions - branch disappeared
                EndOfBranch=1; break; 
         else
               if ~(all(Box(:,1)<=Sol) && all(Sol<=Box(:,2)))
                  OutOfBox=1; break;
               end
               Branch=[Branch,Sol];
               ex0=ex;ex1=Sol;
               ddisp(dflag,['Point=[ ' num2str(ex1',' %0.5g ') ' ]']);
         end
end

Sol.Branch=Branch;
Sol.EndPoints=EndPoints;

Stat.BranchingPoint=BranchingPoint;
Stat.EndOfBranch=EndOfBranch;
Stat.OutOfBox=OutOfBox;

function Sol=SphereStep(SolverFun,SolverOpt,SystemData,Sphere,ToBeRemoved,ScaleCube,SamePoints)

nex=size(Sphere.Center,1);

Limits(:,1)=Sphere.Center-ScaleCube*Sphere.Radius*ones(nex,1);
Limits(:,2)=Sphere.Center+ScaleCube*Sphere.Radius*ones(nex,1);

Sol=feval(SolverFun,@ExtendedSystem,SystemData,Sphere,Limits,SolverOpt);
% Remove previously computed point from new solution set
if ~isempty(ToBeRemoved)
   for i=1:size(Sol,2)
       if max(abs(Sol(:,i)-ToBeRemoved))<=SamePoints
           Sol(:,i)=[];
           break; 
       end
   end
end

function dx=ExtendedSystem(ex,k,Data,Sphere)

if k<size(ex,1)
   dx=feval(Data.Func,ex(1:end-1),k,Data.FuncData,ex(end));
else
   dex=ex-Sphere.Center;
   dx=dex'*dex-Sphere.Radius^2;
end

function ddisp(Flag,varargin)

if Flag
   if nargin<3
      disp(varargin{1});
   else
      disp(sprintf(varargin{1},varargin{2}));
   end
end

function y=nrm(x)

y=x./sqrt(x'*x);
