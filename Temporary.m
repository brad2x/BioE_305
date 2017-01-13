urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
         'part_composition_setup.m');
rehash;
part_composition_setup('v5');

k_dR1 = .0277;
k_dR2 = .0277;
k_dR3 = .0277;
k_dR4 = .0277;

k_R3 = 500*k_dR3;
k_R4 = 500*k_dR4;

K_R1 = 350;
K_R2 = 350;

n_R1 = 3;
n_R2 = 3;

comparator=BioSystem();

comparator.AddConstant('k_dR3', k_dR3);
comparator.AddConstant('k_dR4', k_dR4);

comparator.AddConstant('k_R3', k_R3);
comparator.AddConstant('k_R4', k_R4);

comparator.AddConstant('K_R1',K_R1);
comparator.AddConstant('K_R2',K_R2);

comparator.AddConstant('n_R1', n_R1);
comparator.AddConstant('n_R2', n_R2);

dR1dt = comparator.AddCompositor('R1', 0);
dR2dt = comparator.AddCompositor('R2', 0);
dR3dt = comparator.AddCompositor('R3', 0);
dR4dt = comparator.AddCompositor('R4', 0);
    
InputR1 = @(t) 500*sin((pi/100)*t) + 500;
DerivInputR1 = @(t) (500*pi/100)*cos((pi/100)*t);

InputR2 = @(t) 100;
DerivInputR2 = @(t) 0;

comparator.AddPart(Part('Input 1', dR1dt,...
    Rate('DerivInputR1(t)')));

comparator.AddPart(Part('Input 2', dR2dt,...
    Rate('DerivInputR2(t)')));

comparator.AddPart(Part('NOT1', dR4dt,...
    Rate('k_R4*(((K_R1)^(n_R1))/(((K_R1)^(n_R1))+((InputR1(t))^(n_R1)))) - k_dR4*R4')));

comparator.AddPart(Part('NOT2', dR3dt,...
    Rate('k_R3*(((K_R2)^(n_R2))/(((K_R2)^(n_R2))+((InputR2(t))^(n_R2)))) - k_dR3*R3')));

[T,Y] = comparator.run([0,500]);

% [T, Y] = comparator.run_pulses([...
%    Pulse(0,'R1',100), ...
%    Pulse(200, 'R1', 1000), ... 
%    Pulse(400,'R2',1000), ...
%    Pulse(600,'R1',100), ...
%    Pulse(800,'R2',100)...
%    Pulse(1000,'R1',0)
%    ]);

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

% figure()
% plot(T,Y(:,comparator.CompositorIndex('R1')),'m',T,Y(:,comparator.CompositorIndex('R2')),'c',T,Y(:,comparator.CompositorIndex('GFP')),'g',T,upper_threshold,'r--','LineWidth',2)
% xlabel('Minutes'); ylabel('Molecules/Cell'); title('Magnitude Comparator: GFP v Inputs')
% legend('R1','R2','GFP')
% 
% figure()
% plot(T,Y(:,comparator.CompositorIndex('BFP')),'b',T,Y(:,comparator.CompositorIndex('RFP')),'r',T,Y(:,comparator.CompositorIndex('GFP')),'g',T,upper_threshold,'k--','LineWidth',2)
% xlabel('Minutes'); ylabel('Molecules/Cell'); title('Magnitude Comparator')
% legend('BFP','RFP','GFP','Measurement Threshold')