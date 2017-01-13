urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
         'part_composition_setup.m');
rehash;
part_composition_setup('v5');

k_dR1 = log(2)/25;
k_dR2 = log(2)/25;
k_dR3 = log(2)/5;
k_dR4 = log(2)/5;
k_dR5 = log(2)/5;
k_dBFP = log(2)/25;
k_dRFP = log(2)/25;
k_dGFP = log(2)/25;

k_R1 = 300*k_dR1;
k_R2 = 300*k_dR2;
k_R3 = 500*k_dR3;
k_R4 = 500*k_dR4;
k_R5 = 500*k_dR5;
k_BFP = 300*k_dBFP;
k_RFP = 300*k_dRFP;
k_GFP = 300*k_dGFP;

K_R1 = 350;
K_R2 = 350;
K_R1R3 = K_R1;
K_R3R1 = 150;
K_R2R4 = K_R2;
K_R4R2 = 150;
K_R5 = 250;

n_R1 = 3;
n_R2 = 3;
n_R3 = 3;
n_R4 = 3;
n_R5 = 3;

comparator=BioSystem();
comparator.AddConstant('k_R3', k_R3);
comparator.AddConstant('k_R4', k_R4);
comparator.AddConstant('k_R5', k_R5);
comparator.AddConstant('k_BFP', k_BFP);
comparator.AddConstant('k_RFP', k_RFP);
comparator.AddConstant('k_GFP', k_GFP);

comparator.AddConstant('k_dR1', k_dR1);
comparator.AddConstant('k_dR2', k_dR2);
comparator.AddConstant('k_dR3', k_dR3);
comparator.AddConstant('k_dR4', k_dR4);
comparator.AddConstant('k_dR5', k_dR5);
comparator.AddConstant('k_dBFP', k_dBFP);
comparator.AddConstant('k_dRFP', k_dRFP);
comparator.AddConstant('k_dGFP', k_dGFP);

comparator.AddConstant('K_R1',K_R1);
comparator.AddConstant('K_R2',K_R2);
comparator.AddConstant('K_R1R3',K_R1R3);
comparator.AddConstant('K_R3R1',K_R3R1);
comparator.AddConstant('K_R2R4',K_R2R4);
comparator.AddConstant('K_R4R2',K_R4R2);
comparator.AddConstant('K_R5',K_R5);

comparator.AddConstant('n_R1', n_R1);
comparator.AddConstant('n_R2', n_R2);
comparator.AddConstant('n_R3', n_R3);
comparator.AddConstant('n_R4', n_R4);
comparator.AddConstant('n_R5', n_R5);

dR1dt = comparator.AddCompositor('R1', 0);
dR2dt = comparator.AddCompositor('R2', 0);
dR3dt = comparator.AddCompositor('R3', 450);
dR4dt = comparator.AddCompositor('R4', 450);
dR5dt = comparator.AddCompositor('R5', 0);
dBFPdt = comparator.AddCompositor('BFP', 0);
dRFPdt = comparator.AddCompositor('RFP', 0);
dGFPdt = comparator.AddCompositor('GFP', 0);
    
%comparator.AddPart(Part('R1', [dR1dt],[Rate('-k_dR1*R1')]));

%comparator.AddPart(Part('R2', [dR2dt],[Rate('-k_dR2*R2')]));

comparator.AddPart(Part('NOT1', dR4dt,...
    Rate('k_R4*( ((K_R1)^(n_R1))/( ((K_R1)^(n_R1))+((R1)^(n_R1)) ) ) - k_dR4*R4')));

comparator.AddPart(Part('NOT2', dR3dt,...
    Rate('k_R3*( ((K_R2)^(n_R2))/( ((K_R2)^(n_R2))+((R2)^(n_R2)) ) ) - k_dR3*R3')));

comparator.AddPart(Part('NOR1', [dR5dt, dBFPdt],...
    [Rate('k_R5*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) )'), ...
    Rate('k_BFP*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) ) - k_dBFP*BFP')]));

comparator.AddPart(Part('NOR2', [dR5dt, dRFPdt],...
    [Rate('k_R5*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dR5*R5'), ...
    Rate('k_RFP*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dRFP*RFP')]));

comparator.AddPart(Part('NOT3', dGFPdt,...
    Rate('k_GFP*(((K_R5^n_R5)/((K_R5^n_R5)+(R5^n_R5)))^2) - k_dGFP*GFP')));

%[T,Y] = comparator.run([0,500], odeset('RelTol', 1e-5));

[T, Y] = comparator.run_pulses([...
   Pulse(0,'R1',0), ...
   Pulse(200, 'R1', 1000), ... 
   Pulse(400,'R2',1000), ...
   Pulse(600,'R1',0), ...
   Pulse(800,'R2',0)...
   Pulse(1000,'R1',0)
   ]);

lower_threshold = 25*ones(length(T));
upper_threshold = 275*ones(length(T));

%PLOT ALL GATES REAL-TIME
figure()
subplot(3,2,1)
plot(T,Y(:,1),'m',T,Y(:,4),'y','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('Not1 : R4=NOT(R1)'); legend('R1','R4')

subplot(3,2,2)
plot(T,Y(:,2),'c',T,Y(:,3),'y','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('Not2 : R3=NOT(R2)'); legend('R2','R3')

subplot(3,2,3)
plot(T,Y(:,1),'m',T,Y(:,3),'y',T,Y(:,6),'b','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOR1: R5=NOR(R1,R3)'); legend('R1','R3','BFP')

subplot(3,2,4)
plot(T,Y(:,2),'c',T,Y(:,4),'y',T,Y(:,7),'r','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOR2: RFP=NOR(R2,R4)'); legend('R2','R4','RFP')

subplot(3,2,5)
plot(T,Y(:,5),'k',T,Y(:,8),'g',T,upper_threshold,'r--','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOT3: GFP=NOT(R5)'); legend('R5','GFP')

subplot(3,2,6)
plot(T,Y(:,1),'m',T,Y(:,2),'c',T,Y(:,6),'b',T,Y(:,7),'r',T,Y(:,8),'g',T,upper_threshold,'k--')
xlabel('Minutes'); ylabel('Molecules/Cell');
title('Comparator: R:R1<R2 B:R1>R2 G:R1=R2'); 
legend('R1','R2','BFP','RFP','GFP')

figure()
plot(T,Y(:,comparator.CompositorIndex('R1')),'m',T,Y(:,comparator.CompositorIndex('R2')),'c',T,Y(:,comparator.CompositorIndex('GFP')),'g',T,upper_threshold,'r--','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell'); title('Magnitude Comparator: GFP v Inputs')
legend('R1','R2','GFP')

figure()
plot(T,Y(:,comparator.CompositorIndex('BFP')),'b',T,Y(:,comparator.CompositorIndex('RFP')),'r',T,Y(:,comparator.CompositorIndex('GFP')),'g',T,upper_threshold,'k--','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell'); title('Magnitude Comparator')
legend('BFP','RFP','GFP','Measurement Threshold')