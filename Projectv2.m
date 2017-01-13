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

comparator = BioSystem();

comparator.AddConstant('k_pR1', k_pR1);
comparator.AddConstant('k_pR2', k_pR2);
comparator.AddConstant('k_pR3', k_pR3);
comparator.AddConstant('k_pR4', k_pR4);
comparator.AddConstant('k_pBFP', k_pBFP);
comparator.AddConstant('k_pRFP', k_pRFP);
comparator.AddConstant('k_pGFP', k_pGFP);

comparator.AddConstant('k_degA', k_degA);
comparator.AddConstant('k_degB', k_degB);
comparator.AddConstant('k_degR1', k_degR1);
comparator.AddConstant('k_degR2', k_degR2);
comparator.AddConstant('k_degR3', k_degR3);
comparator.AddConstant('k_degR4', k_degR4);
comparator.AddConstant('k_degBFP', k_degBFP);
comparator.AddConstant('k_degRFP', k_degRFP);
comparator.AddConstant('k_degGFP', k_degGFP);

comparator.AddConstant('K_A',K_A);
comparator.AddConstant('K_B',K_B);
comparator.AddConstant('K_R1',K_R1);
comparator.AddConstant('K_R2',K_R2);
comparator.AddConstant('K_R3',K_R3);
comparator.AddConstant('K_R4',K_R4);

comparator.AddConstant('n_A', n_A);
comparator.AddConstant('n_B', n_B);
comparator.AddConstant('n_R1', n_R1);
comparator.AddConstant('n_R2', n_R2);
comparator.AddConstant('n_R3', n_R3);
comparator.AddConstant('n_R4', n_R4);

dAdt = comparator.AddCompositor('A', 0);
dBdt = comparator.AddCompositor('B', 0);
dR1dt = comparator.AddCompositor('R1', 10);
dR2dt = comparator.AddCompositor('R2', 10);
dR3dt = comparator.AddCompositor('R3', 10);
dR4dt = comparator.AddCompositor('R4', 10);
dBFPdt = comparator.AddCompositor('BFP', 0);
dRFPdt = comparator.AddCompositor('RFP', 0);
dGFPdt = comparator.AddCompositor('GFP', 0);

comparator.AddPart(Part('A->R1', [dAdt dR1dt],...
    [Rate('-k_degA*A'), Rate('k_pR1*( (K_A^n_A)/( (K_A^n_A)+(A^n_A) ) ) - k_degR1*R1') ]));

comparator.AddPart(Part('B->R2', [dBdt dR2dt],...
    [Rate('-k_degB*B'), Rate('k_pR2*( (K_B^n_B)/( (K_B^n_B)+(B^n_B) ) ) - k_degR2*R2') ]));

comparator.AddPart(Part('A_+ R2_-> R3 + BFP', [dAdt dR2dt dR3dt dBFPdt],...
    [Rate('0'), Rate('0'),...
    Rate('k_pR3*(((K_A^n_A)/((K_A^n_A)+(A^n_A)))*((K_R2^n_R2)/((K_R2^n_R2)+(R2^n_R2)))) - k_degR3*R3'),...
    Rate('k_pBFP*(((K_A^n_A)/((K_A^n_A)+(A^n_A)))*((K_R2^n_R2)/((K_R2^n_R2)+(R2^n_R2)))) - k_degBFP*BFP')]));

comparator.AddPart(Part('B_+R1_-> R4 + RFP', [dBdt dR1dt dR4dt dRFPdt],...
    [Rate('0'), Rate('0'),...
    Rate('k_pR4*(((K_B^n_B)/((K_B^n_B)+(B^n_B)))*((K_R1^n_R1)/((K_R1^n_R1)+(R1^n_R1)))) - k_degR4*R4'),...
    Rate('k_pRFP*(((K_B^n_B)/((K_B^n_B)+(B^n_B)))*((K_R1^n_R1)/((K_R1^n_R1)+(R1^n_R1)))) - k_degRFP*RFP')]));

comparator.AddPart(Part('R3_+R4_-> GFP', [dR3dt dR4dt dGFPdt],...
    [Rate('0'), Rate('0'),...
    Rate('k_pGFP*(((K_R3^n_R3)/((K_R3^n_R3)+(R3^n_R3)))*((K_R4^n_R4)/((K_R4^n_R4)+(R4^n_R4)))) - k_degGFP*GFP')]));


%[T,Y] = comparator.run([0,100]);

[T, Y] = comparator.run_pulses([...
    Pulse(0, 'A', 10), Pulse(.0001,'B',10),... 
    Pulse(1000, 'A', 1000), ...
    Pulse(2000, 'A', 0), ... 
    Pulse(3000, 'B',1000), ...
    Pulse(4000,'B',0), ...
    Pulse(5000,'A',1000), Pulse(5000.0001,'B',1000),...
    Pulse(6000,'A',0), Pulse(6000.0001,'B',0),...
    Pulse(7000,'A',0)
]);

threshold = ones(length(T))*5000;

figure()
subplot(2,3,1)
plot(T,Y(:,comparator.CompositorIndex('A')),T,Y(:,comparator.CompositorIndex('R1')))
xlabel('time'); ylabel('concentration'); legend('A','R1'); title('R1 = NOT(A)');

subplot(2,3,2)
plot(T,Y(:,comparator.CompositorIndex('B')),T,Y(:,comparator.CompositorIndex('R2')))
xlabel('time'); ylabel('concentration'); legend('B','R2'); title('R2 = NOT(B)');

subplot(2,3,3)
plot(T,Y(:,comparator.CompositorIndex('A')),T,Y(:,comparator.CompositorIndex('R2')),...
    T,Y(:,comparator.CompositorIndex('R3')),T,Y(:,comparator.CompositorIndex('BFP')))
xlabel('time'); ylabel('concentration'); legend('A','R2','R3','BFP');
title('R3 = BFP = NOR(A,R2)');

subplot(2,3,4)
plot(T,Y(:,comparator.CompositorIndex('B')),T,Y(:,comparator.CompositorIndex('R1')),...
    T,Y(:,comparator.CompositorIndex('R4')),T,Y(:,comparator.CompositorIndex('RFP')))
xlabel('time'); ylabel('concentration'); legend('B','R1','R4','RFP');
title('R4 = RFP = NOR(B,R1)');

subplot(2,3,5)
plot(T,Y(:,comparator.CompositorIndex('R3')),T,Y(:,comparator.CompositorIndex('R4')),...
    T,Y(:,comparator.CompositorIndex('GFP')))
xlabel('time'); ylabel('concentration'); legend('R3','R4','GFP');
title('GFP = NOR(R3,R4)');

figure()
subplot(2,2,1)
plot(T,Y(:,comparator.CompositorIndex('A')),T,Y(:,comparator.CompositorIndex('B')),...
    T,Y(:,comparator.CompositorIndex('RFP')),T,threshold)
xlabel('time'); ylabel('concentration'); legend('A','B','RFP','Threshold');
title('RFP = A>B');

subplot(2,2,2)
plot(T,Y(:,comparator.CompositorIndex('A')),T,Y(:,comparator.CompositorIndex('B')),...
    T,Y(:,comparator.CompositorIndex('BFP')),T,threshold)
xlabel('time'); ylabel('concentration'); legend('A','B','GFP','Threshold');
title('BFP = A<B)');

subplot(2,2,3)
plot(T,Y(:,comparator.CompositorIndex('A')),T,Y(:,comparator.CompositorIndex('B')),...
    T,Y(:,comparator.CompositorIndex('GFP')),T,threshold)
xlabel('time'); ylabel('concentration'); legend('A','B','GFP','Threshold');
title('GFP = A=B');

subplot(2,2,4)
plot(T,Y(:,comparator.CompositorIndex('A')),T,Y(:,comparator.CompositorIndex('B')),...
    T,Y(:,comparator.CompositorIndex('R3')),T,Y(:,comparator.CompositorIndex('R4')))
xlabel('time'); ylabel('concentration'); legend('A','B','R3','R4');
title('Intermediates');

%GATE PLOTTING

output = ones(1001);

for A_Val = 0:1:1000
    not1=BioSystem();
    not1.AddConstant('k_pR1', k_pR1);
    not1.AddConstant('k_degA', 0);
    not1.AddConstant('k_degR1', k_degR1);
    not1.AddConstant('K_A',K_A);
    not1.AddConstant('n_A', n_A);
    
    dAdt = not1.AddCompositor('A', A_Val);
    dR1dt = not1.AddCompositor('R1', 10);
    
    not1.AddPart(Part('A->R1', [dAdt dR1dt],...
    [Rate('-k_degA*A'), Rate('k_pR1*( (K_A^n_A)/( (K_A^n_A)+(A^n_A) ) ) - k_degR1*R1') ]));
    
    [T,Y] = comparator.run([0,500]);
    output(A_Val)=Y(end,not1.CompositorIndex('R1'));
end

figure()
plot(A_Val,output) 
xlabel('[A]'); ylabel('R1'); title('Not1 Transfer Fxn')
    