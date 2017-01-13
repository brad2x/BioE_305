A = 500;
w = (2*pi)/150;
C = 500;

%osc = A*sin(wt)+A
%osc' = A*w*cos(wt)

sys = BioSystem();

sys.AddConstant('A', A);
sys.AddConstant('w', w);
sys.AddConstant('C', C);

dtimedt = sys.AddCompositor('T1',0);
dR1dt = sys.AddCompositor('R1',A);

sys.AddPart(Part('time', dtimedt, Rate('1')));

sys.AddPart(Part('R1', dR1dt, Rate('A*w*cos(w*T1)')));

[T,Y] = sys.run([0,1000]);

plot(T,Y(:,2));