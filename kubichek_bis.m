function kubichek_bis(cp)

par.cp=cp;
par.A=2;par.B=6;par.R0=0.1;

% GridMachine parameters
opt.Limits=[zeros(8,1),10*ones(8,1)];
opt.StackSize=500000;
opt.InitialGrid=2;
opt.ZeroCube=1e-3;
opt.Progress=1;

[Sol,Stat]=GridMachineBis(@kubicheck,par,opt);

function f=kubicheck(X,k,par)

A=par.A;B=par.B;R0=par.R0;cp=par.cp;

NPoints=size(X,2);
f=zeros(1,NPoints);

for i=1:NPoints
    
    X1=X(1,i);  Y1=X(2,i);  X2=X(3,i);  Y2=X(4,i);
    X3=X(5,i);  Y3=X(6,i);  X4=X(7,i);  Y4=X(8,i);
    D1=10^cp; 	D2=D1/R0;

    if k==1
       f(i)=A-(B+1)*X1+X1*X1*Y1+D1*(X2-X1);
    elseif k==2
       f(i)=B*X1-X1*X1*Y1+D2*(Y2-Y1);
    elseif k==3
       f(i)=A-(B+1)*X2+X2*X2*Y2+D1*(X1-2*X2+X3);
    elseif k==4
       f(i)=B*X2-X2*X2*Y2+D2*(Y1-2*Y2+Y3);
    elseif k==5
       f(i)=A-(B+1)*X3+X3*X3*Y3+D1*(X2-2*X3+X4);
    elseif k==6
       f(i)=B*X3-X3*X3*Y3+D2*(Y2-2*Y3+Y4);
    elseif k==7
       f(i)=A-(B+1)*X4+X4*X4*Y4+D1*(X3-X4);
    elseif k==8
       f(i)=B*X4-X4*X4*Y4+D2*(Y3-Y4);
    end
end
