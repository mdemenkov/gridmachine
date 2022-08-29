function out=test_C


out=ex_kubichek;

function out=ex_sphere

par.Limits=[-3 3;-3 3;-3 3];
par.DeriSteps=[0.05;0.05;0.05];
par.MaxTime=50;
par.NewtonZero=1e-5;
par.JacStep=1;
par.Density=5;
par.Radius=4;
par.SolSize=10;
par.SamePoints=1e-2;

out=GridMachineC(par);

function out=ex_kubichek

par.Limits=[zeros(8,1),10*ones(8,1)];
par.Density=3;
par.DeriSteps=0.02*ones(8,1);
par.MaxTime=50;
par.NewtonZero=1e-5;
par.JacStep=1;
par.Radius=0;
par.SolSize=20;
par.SamePoints=1e-2;

out=GridMachineC(par);
