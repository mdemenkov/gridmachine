function Sol=test_5order(u)

%	The peculiarities of the aircraft spatial motion in the case of
%	---------------------------------------------------------------
%	rapid rolling.
%	--------------
%
%	Ref.: M.G.Goman "Differential Method of Continuation of the
%	Solutions of the Systems of Finite Nonlinear Equations Depen-
%	ding on the Parameter", Uchenye Zapiski TsAGI, vol. XVII, no.5,
%	1986. (NB There are some mistakes in the paper concerning the
%	formulation of the problem)

% Task parameters
A=zeros(5,5);
A(1,1)=-.7;A(2,2)=-.17;A(3,1)=-4./2.;A(3,3)=-.8/2;
A(4,2)=-3.*2.;A(4,4)=-.6;A(4,5)=.055;A(5,2)=-15.5;
A(5,4)=-.7;A(5,5)=-2.35;
par.A=A;
B=zeros(5,3);
B(1,1)=-.002;B(2,3)=-.0004;B(3,1)=-.156;
B(4,2)=.0035;B(4,3)=-.05;
B(5,2)=-.5;B(5,3)=-.07;
par.B=B;

par.C1=-.9;par.C2=.8;

par.DB=u(1);% stabilizer deflection, deg.
par.DE=u(2);% aileron    deflection, deg.
par.DN=u(2);% rudder     deflection, deg.

% GridMachine parameters
opt_gm.Limits=[-1.57 1.57;-1 1;-.5 .5;-1 1;-5 5];
opt_gm.Density=3;
opt_gm.DeriSteps=0.05*ones(5,1);
opt_gm.Progress=1;

[Sol,Stat]=GridMachine(@article_example,par,opt_gm)


function f=article_example(X,k,par)

% State variables:
% x(1) \alpha angle of attack, radians
% x(2) \beta sideslip angle,  radians
% x(3) w_z pitch rate, 1/sec
% x(4) w_y yaw   rate, 1/sec
% x(5) w_x roll  rate, 1/sec
        
A=par.A;B=par.B;C1=par.C1;C2=par.C2;
DB=par.DB;DE=par.DE;DN=par.DN;

N=size(X,2);
f=zeros(1,N);
for i=1:N
    x=X(:,i);
    if k==1
       f(i)=A(1,1)*x(1)+x(3)+B(1,1)*DB-x(2)*x(5);
    elseif k==2
       f(i)=A(2,2)*x(2)+x(4)+B(2,3)*DN+x(1)*x(5);
    elseif k==3
       f(i)=A(3,1)*x(1)+A(3,3)*x(3)+B(3,1)*DB+C1*x(4)*x(5);
    elseif k==4
       f(i)=A(4,2)*x(2)+A(4,4)*x(4)+A(4,5)*x(5)+B(4,2)*DE+B(4,3)*DN+C2*x(5)*x(3);
    elseif k==5
       f(i)=A(5,2)*x(2)+A(5,4)*x(4)+A(5,5)*x(5)+B(5,2)*DE+B(5,3)*DN;
    end
end
