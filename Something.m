urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
         'part_composition_setup.m');
rehash;
part_composition_setup('v5');

k_pR1 = 1000;
k_pR2 = 1000;
k_pR3 = 1000;
k_pR4 = 1000;
k_pBFP = 100;
k_pRFP = 100;
k_pGFP = 100;

k_degA = .001;
k_degB = .001;
k_degR1 = .1;
k_degR2 = .1;
k_degR3 = .1;
k_degR4 = .1;
k_degBFP = .01;
k_degRFP = .01;
k_degGFP = .01;

K_A = 10;
K_B = 10;
K_R1 = 10;
K_R2 = 10;
K_R3 = 20;
K_R4 = 20;

n_A = 2;
n_B = 2;
n_R1 = 2;
n_R2 = 2;
n_R3 = 2;
n_R4 = 2;

range = 100;
output = zeros(length(0:1:range));
input = zeros(length(0:1:range));
n = 0;

for A_Val = 0:1:range
    n = n+1;
    not1=BioSystem();
    not1.AddConstant('k_pR1', k_pR1);
    not1.AddConstant('k_degA', 0);
    not1.AddConstant('k_degR1', k_degR1);
    not1.AddConstant('K_A',K_A);
    not1.AddConstant('n_A', n_A);
    
    dAdt = not1.AddCompositor('A', A_Val);
    dR1dt = not1.AddCompositor('R1', 0);
    
    not1.AddPart(Part('A->R1', [dAdt dR1dt],...
    [Rate('-k_degA*A'), Rate('k_pR1*( (K_A^n_A)/( (K_A^n_A)+(A^n_A) ) ) - k_degR1*R1') ]));
    
    [T,Y] = not1.run([0,500]);
    input(n) = A_Val;
    output(n)=Y(end,not1.CompositorIndex('R1'));
end

figure()
plot(input,output) 
xlabel('[A]'); ylabel('R1'); title('Not1 Transfer Fxn')