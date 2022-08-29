function test_wing_bis

   % System parameters
   coef.dBeta=0.2;
   coef.k=0.008;coef.theta_0=0.5;coef.C_lb0=-0.4;coef.C_lp0=-0.1;
   coef.tau_1=10;

   Limits=[-2 2;-2 2;-2 2];
   opt.InitialGrid=10;
   opt.StackSize=1000;
   opt.Progress=0;   
    
   % Compute initial equilibriums
   cp0=4;
   Sol=BisectionMethod(@wing_rhs,coef,cp0,Limits,opt);
   x0=Sol(:,1);
   
   % Start continuation at [x0;cp0]
   BisOpt.InitialGrid=2;
   BisOpt.StackSize=100000;
   BisOpt.ZeroCube=0.00001;
   BisOpt.SamePoints=0.001;
   BisOpt.Progress=1;
   
   SphereOpt.SolverFun=@BisectionMethod;
   SphereOpt.SolverOpt=BisOpt;
   SphereOpt.StackSize=10;
   SphereOpt.Box=[-2 2;-2 2;-2 2;-4 4];
   SphereOpt.MaxRad=0.2;
   SphereOpt.MinRad=0.02;
   SphereOpt.xcp0=[x0;cp0];
   SphereOpt.ScaleCube=2;
   SphereOpt.ScaleSphere=1/2;
   SphereOpt.SamePoints=0.001;
   SphereOpt.Progress=1;
 
   Branches=SphereMethod(@wing_rhs,coef,SphereOpt);

   % Plot results
   figure;hold on;xlabel('k_2');ylabel('\beta');
   axis([-4 4 -2 2]);
   for i=1:length(Branches)
       Branch=Branches{i};
       plot(Branch(4,:),Branch(1,:),'b');
   end

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
