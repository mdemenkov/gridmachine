function test_wing65

%-Flags------------
ComputeDirect=0;
ComputeBranches=1;
DisplayResults=1;
%------------------

coef.dBeta=0.2;
coef.BetaMin=-4;coef.BetaMax=4;
coef.k=0.008;coef.theta_0=0.5;coef.C_lb0=-0.4;coef.C_lp0=-0.1;
coef.tau_1=10;

k2_range=-4:0.2:4;

Limits=[-2 2;-2 2;-2 2];
opt.Density=10;
opt.DeriSteps=[0.01;0.01;0.01];
opt.Progress=0;   

if ComputeDirect
    j=0;    
    for k2=k2_range
        disp(['Computing for k_2=' num2str(k2) ' ...']);
        States=NewtonGrid(@wing_rhs,coef,k2,Limits,opt);
        for i=1:size(States,2)
            j=j+1;
            X(j)=k2;Y(j)=States(1,i);COL(j,:)='b.';
        end
    end
    figure;hold on;xlabel('k_2');ylabel('\beta');
    for i=1:j
        plot(X(i),Y(i),char(COL(i,:)));
    end    
end

if ComputeBranches
    
   cp0=4;
   States=NewtonGrid(@wing_rhs,coef,cp0,Limits,opt);
   x0=States(:,1);
   
   NGOpt.Density=2;
   NGOpt.DeriSteps=0.001*ones(4,1);
   NGOpt.Progress=0;
   
   BMOpt.SolverFun=@NewtonGrid;
   BMOpt.SolverOpt=NGOpt;
   BMOpt.StackSize=10;
   BMOpt.Box=[-2 2;-2 2;-2 2;-4 4];
   BMOpt.MaxRad=0.2;
   BMOpt.MinRad=0.02;
   BMOpt.xcp0=[x0;cp0];
   BMOpt.ScaleCube=2;
   BMOpt.ScaleSphere=1/2;
   BMOpt.SamePoints=0.001;
   BMOpt.Progress=1;
tic   
   Branches=SphereMethod(@wing_rhs,coef,BMOpt);
toc
end

if DisplayResults
   figure;hold on;xlabel('k_2');ylabel('\beta');
   axis([-4 4 -2 2]);
   for i=1:length(Branches)
       Branch=Branches{i};
       plot(Branch(4,:),Branch(1,:),'b');
   end
end


% RHS ------------------------------

function dx=wing_rhs(x,k,coef,k2)
% Model of a free-to-roll delta wing mounted on a fixed  sting.
% Reference: M.G. Goman, A.N. Krabrov, A.V. Khramtsovsky, 
% "Chaotic dynamics in a simple aeromechanical system", 
% in J.M. Blackledge, A.K. Evans and M.J. Turner (Eds.) 
% "Fractal geometry: mathematical methods, algorithms, applications",
% Horwood Publishing, Chichester, 2002, pp. 1-15.

% k2 is the continuation parameter
beta=x(1);p=x(2);C_lv=x(3);
% deform=0.01;
deform=0;
switch k
       case 1 % dbeta
            dx=p*sin(coef.theta_0);
       case 2 % dp
            dx=coef.k*(coef.C_lb0*beta+coef.C_lp0*p+C_lv);
       case 3 % dC_lv
            dx=(-C_lv+deform)/coef.tau_1+k2*coef.C_lb0*beta/(1+beta^2)/coef.tau_1;
end
