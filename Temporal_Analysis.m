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

parameter_start = .1;
parameter_range = 25;
parameter_step = 1;

threshold_time_BFP = zeros(length(parameter_start:parameter_step:parameter_range),1);
parameter_input_BFP = zeros(length(parameter_start:parameter_step:parameter_range),1);

m = 0;
for parameter = parameter_start:parameter_step:parameter_range
    m = m+1;

    R3deg=BioSystem();
    R3deg.AddConstant('k_R3', 500*log(2)/parameter);
    R3deg.AddConstant('k_R4', k_R4);
    R3deg.AddConstant('k_R5', k_R5);
    R3deg.AddConstant('k_BFP', k_BFP);
    R3deg.AddConstant('k_RFP', k_RFP);
    R3deg.AddConstant('k_GFP', k_GFP);

    R3deg.AddConstant('k_dR3', log(2)/parameter);
    R3deg.AddConstant('k_dR4', k_dR4);
    R3deg.AddConstant('k_dR5', k_dR5);
    R3deg.AddConstant('k_dBFP', k_dBFP);
    R3deg.AddConstant('k_dRFP', k_dRFP);
    R3deg.AddConstant('k_dGFP', k_dGFP);

    R3deg.AddConstant('K_R1',K_R1);
    R3deg.AddConstant('K_R2',K_R2);
    R3deg.AddConstant('K_R1R3',K_R1R3);
    R3deg.AddConstant('K_R3R1',K_R3R1);
    R3deg.AddConstant('K_R2R4',K_R2R4);
    R3deg.AddConstant('K_R4R2',K_R4R2);
    R3deg.AddConstant('K_R5',K_R5);

    R3deg.AddConstant('n_R1', n_R1);
    R3deg.AddConstant('n_R2', n_R2);
    R3deg.AddConstant('n_R3', n_R3);
    R3deg.AddConstant('n_R4', n_R4);
    R3deg.AddConstant('n_R5', n_R5);

    dR1dt = R3deg.AddCompositor('R1', 0);
    dR2dt = R3deg.AddCompositor('R2', 0);
    dR3dt = R3deg.AddCompositor('R3', 500);
    dR4dt = R3deg.AddCompositor('R4', 500);
    dR5dt = R3deg.AddCompositor('R5', 0);
    dBFPdt = R3deg.AddCompositor('BFP', 0);
    dRFPdt = R3deg.AddCompositor('RFP', 0);
    dGFPdt = R3deg.AddCompositor('GFP', 300);   

    R3deg.AddPart(Part('NOT1', dR4dt,...
        Rate('k_R4*( ((K_R1)^(n_R1))/( ((K_R1)^(n_R1))+((R1)^(n_R1)) ) ) - k_dR4*R4')));

    R3deg.AddPart(Part('NOT2', dR3dt,...
        Rate('k_R3*( ((K_R2)^(n_R2))/( ((K_R2)^(n_R2))+((R2)^(n_R2)) ) ) - k_dR3*R3')));

    R3deg.AddPart(Part('NOR1', [dR5dt, dBFPdt],...
        [Rate('k_R5*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) )'), ...
        Rate('k_BFP*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) ) - k_dBFP*BFP')]));

    R3deg.AddPart(Part('NOR2', [dR5dt, dRFPdt],...
        [Rate('k_R5*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dR5*R5'), ...
        Rate('k_RFP*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dRFP*RFP')]));

    R3deg.AddPart(Part('NOT3', dGFPdt,...
        Rate('k_GFP*(((K_R5^n_R5)/((K_R5^n_R5)+(R5^n_R5)))^2) - k_dGFP*GFP')));

    [T_BFP,Y_BFP] = R3deg.run_pulses([...
    Pulse(0,'R2',0), ...
    Pulse(100,'R2',1000), ...
    Pulse(800,'R2',0) ...
    ]);
    
    upper_threshold_BFP = 275*ones(length(T_BFP));
    parameter_input_BFP(m) = parameter;
    
    min_above_threshold = find(Y_BFP(:,R3deg.CompositorIndex('BFP')) >= upper_threshold_BFP,1);
    if( isempty(min_above_threshold) )
        min_above_threshold = length(T_BFP);
    end
    threshold_time_BFP(m) = T_BFP(min_above_threshold) - 100;
end

threshold_time_RFP = zeros(length(parameter_start:parameter_step:parameter_range),1);
parameter_input_RFP = zeros(length(parameter_start:parameter_step:parameter_range),1);

m = 0;
for parameter = parameter_start:parameter_step:parameter_range
    m = m+1;

    R4deg=BioSystem();
    R4deg.AddConstant('k_R3', k_R3);
    R4deg.AddConstant('k_R4', 500*log(2)/parameter);
    R4deg.AddConstant('k_R5', k_R5);
    R4deg.AddConstant('k_BFP', k_BFP);
    R4deg.AddConstant('k_RFP', k_RFP);
    R4deg.AddConstant('k_GFP', k_GFP);

    R4deg.AddConstant('k_dR3', k_dR3);
    R4deg.AddConstant('k_dR4', log(2)/parameter);
    R4deg.AddConstant('k_dR5', k_dR5);
    R4deg.AddConstant('k_dBFP', k_dBFP);
    R4deg.AddConstant('k_dRFP', k_dRFP);
    R4deg.AddConstant('k_dGFP', k_dGFP);

    R4deg.AddConstant('K_R1',K_R1);
    R4deg.AddConstant('K_R2',K_R2);
    R4deg.AddConstant('K_R1R3',K_R1R3);
    R4deg.AddConstant('K_R3R1',K_R3R1);
    R4deg.AddConstant('K_R2R4',K_R2R4);
    R4deg.AddConstant('K_R4R2',K_R4R2);
    R4deg.AddConstant('K_R5',K_R5);

    R4deg.AddConstant('n_R1', n_R1);
    R4deg.AddConstant('n_R2', n_R2);
    R4deg.AddConstant('n_R3', n_R3);
    R4deg.AddConstant('n_R4', n_R4);
    R4deg.AddConstant('n_R5', n_R5);

    dR1dt = R4deg.AddCompositor('R1', 0);
    dR2dt = R4deg.AddCompositor('R2', 0);
    dR3dt = R4deg.AddCompositor('R3', 500);
    dR4dt = R4deg.AddCompositor('R4', 500);
    dR5dt = R4deg.AddCompositor('R5', 0);
    dBFPdt = R4deg.AddCompositor('BFP', 0);
    dRFPdt = R4deg.AddCompositor('RFP', 0);
    dGFPdt = R4deg.AddCompositor('GFP', 300);   

    R4deg.AddPart(Part('NOT1', dR4dt,...
        Rate('k_R4*( ((K_R1)^(n_R1))/( ((K_R1)^(n_R1))+((R1)^(n_R1)) ) ) - k_dR4*R4')));

    R4deg.AddPart(Part('NOT2', dR3dt,...
        Rate('k_R3*( ((K_R2)^(n_R2))/( ((K_R2)^(n_R2))+((R2)^(n_R2)) ) ) - k_dR3*R3')));

    R4deg.AddPart(Part('NOR1', [dR5dt, dBFPdt],...
        [Rate('k_R5*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) )'), ...
        Rate('k_BFP*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) ) - k_dBFP*BFP')]));

    R4deg.AddPart(Part('NOR2', [dR5dt, dRFPdt],...
        [Rate('k_R5*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dR5*R5'), ...
        Rate('k_RFP*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dRFP*RFP')]));

    R4deg.AddPart(Part('NOT3', dGFPdt,...
        Rate('k_GFP*(((K_R5^n_R5)/((K_R5^n_R5)+(R5^n_R5)))^2) - k_dGFP*GFP')));

    [T_RFP,Y_RFP] = R4deg.run_pulses([...
    Pulse(0,'R1',0), ...
    Pulse(100,'R1',1000), ...
    Pulse(800,'R1',0) ...
    ]);
    
    upper_threshold_RFP = 275*ones(length(T_RFP));
    parameter_input_RFP(m) = parameter;
    
    min_above_threshold = find(Y_RFP(:,R4deg.CompositorIndex('RFP')) >= upper_threshold_RFP,1);
    if( isempty(min_above_threshold) )
        min_above_threshold = length(T_RFP);
    end
    threshold_time_RFP(m) = T_RFP(min_above_threshold) - 100;
end

threshold_time_GFP = zeros(length(parameter_start:parameter_step:parameter_range),1);
parameter_input_GFP = zeros(length(parameter_start:parameter_step:parameter_range),1);

m = 0;
for parameter = parameter_start:parameter_step:parameter_range
    m = m+1;

    R5deg=BioSystem();
    R5deg.AddConstant('k_R3', k_R3);
    R5deg.AddConstant('k_R4', k_R4);
    R5deg.AddConstant('k_R5', 500*log(2)/parameter);
    R5deg.AddConstant('k_BFP', k_BFP);
    R5deg.AddConstant('k_RFP', k_RFP);
    R5deg.AddConstant('k_GFP', k_GFP);

    R5deg.AddConstant('k_dR3', k_dR3);
    R5deg.AddConstant('k_dR4', k_dR4);
    R5deg.AddConstant('k_dR5', log(2)/parameter);
    R5deg.AddConstant('k_dBFP', k_dBFP);
    R5deg.AddConstant('k_dRFP', k_dRFP);
    R5deg.AddConstant('k_dGFP', k_dGFP);

    R5deg.AddConstant('K_R1',K_R1);
    R5deg.AddConstant('K_R2',K_R2);
    R5deg.AddConstant('K_R1R3',K_R1R3);
    R5deg.AddConstant('K_R3R1',K_R3R1);
    R5deg.AddConstant('K_R2R4',K_R2R4);
    R5deg.AddConstant('K_R4R2',K_R4R2);
    R5deg.AddConstant('K_R5',K_R5);

    R5deg.AddConstant('n_R1', n_R1);
    R5deg.AddConstant('n_R2', n_R2);
    R5deg.AddConstant('n_R3', n_R3);
    R5deg.AddConstant('n_R4', n_R4);
    R5deg.AddConstant('n_R5', n_R5);

    dR1dt = R5deg.AddCompositor('R1', 1000);
    dR2dt = R5deg.AddCompositor('R2', 0);
    dR3dt = R5deg.AddCompositor('R3', 500);
    dR4dt = R5deg.AddCompositor('R4', 0);
    dR5dt = R5deg.AddCompositor('R5', 500);
    dBFPdt = R5deg.AddCompositor('BFP', 0);
    dRFPdt = R5deg.AddCompositor('RFP', 300);
    dGFPdt = R5deg.AddCompositor('GFP', 0);   

    R5deg.AddPart(Part('NOT1', dR4dt,...
        Rate('k_R4*( ((K_R1)^(n_R1))/( ((K_R1)^(n_R1))+((R1)^(n_R1)) ) ) - k_dR4*R4')));

    R5deg.AddPart(Part('NOT2', dR3dt,...
        Rate('k_R3*( ((K_R2)^(n_R2))/( ((K_R2)^(n_R2))+((R2)^(n_R2)) ) ) - k_dR3*R3')));

    R5deg.AddPart(Part('NOR1', [dR5dt, dBFPdt],...
        [Rate('k_R5*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) )'), ...
        Rate('k_BFP*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) ) - k_dBFP*BFP')]));

    R5deg.AddPart(Part('NOR2', [dR5dt, dRFPdt],...
        [Rate('k_R5*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dR5*R5'), ...
        Rate('k_RFP*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dRFP*RFP')]));

    R5deg.AddPart(Part('NOT3', dGFPdt,...
        Rate('k_GFP*(((K_R5^n_R5)/((K_R5^n_R5)+(R5^n_R5)))^2) - k_dGFP*GFP')));

    [T_GFP,Y_GFP] = R5deg.run_pulses([...
    Pulse(0,'R2',0), ...
    Pulse(100,'R2',1000), ...
    Pulse(800,'R2',0) ...
    ]);
    
    upper_threshold_GFP = 275*ones(length(T_GFP));
    parameter_input_GFP(m) = parameter;
    
    min_above_threshold = find(Y_GFP(:,R5deg.CompositorIndex('GFP')) >= upper_threshold_GFP,1);
    if( isempty(min_above_threshold) )
        min_above_threshold = length(T_GFP);
    end
    threshold_time_GFP(m) = T_GFP(min_above_threshold) - 100;
end

figure()
subplot(1,3,1)
plot(parameter_input_BFP, threshold_time_BFP,'bo','LineWidth',2)
xlabel('R3 half life'); ylabel('delay'); 
title('BFP Delay');

subplot(1,3,2)
plot(parameter_input_RFP, threshold_time_RFP,'ro','LineWidth',2)
xlabel('R4 half life'); ylabel('delay'); 
title('RFP Delay');

subplot(1,3,3)
plot(parameter_input_GFP, threshold_time_GFP,'go','LineWidth',2)
xlabel('R5 half life'); ylabel('delay'); 
title('GFP delay');

%% Extra Figures
figure()
subplot(3,2,1)
plot(T_BFP,Y_BFP(:,1),'m',T_BFP,Y_BFP(:,4),'y','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOT1'); legend('R1','R4')

subplot(3,2,2)
plot(T_BFP,Y_BFP(:,2),'c',T_BFP,Y_BFP(:,3),'y','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOT2'); legend('R2','R3')

subplot(3,2,3)
plot(T_BFP,Y_BFP(:,1),'m',T_BFP,Y_BFP(:,3),'y',T_BFP,Y_BFP(:,6),'b',T_BFP,upper_threshold_BFP,'k--','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOR1: BFP=NOR(R1,R3)'); legend('R1','R3','BFP')

subplot(3,2,4)
plot(T_BFP,Y_BFP(:,2),'c',T_BFP,Y_BFP(:,4),'y',T_BFP,Y_BFP(:,7),'r',T_BFP,upper_threshold_BFP,'k--','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOR2: RFP=NOR(R2,R4)'); legend('R2','R4','RFP')

subplot(3,2,5)
plot(T_BFP,Y_BFP(:,5),'k',T_BFP,Y_BFP(:,8),'g',T_BFP,upper_threshold_BFP,'r--','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOT3: GFP=NOT(R5)'); legend('R5','GFP')

subplot(3,2,6)
plot(T_BFP,Y_BFP(:,1),'m',T_BFP,Y_BFP(:,2),'c',T_BFP,Y_BFP(:,6),'b',T_BFP,Y_BFP(:,7),'r',T_BFP,Y_BFP(:,8),'g',T_BFP,upper_threshold_BFP,'k--')
xlabel('Minutes'); ylabel('Molecules/Cell');
title('Comparator: R:R1<R2 B:R1>R2 G:R1=R2'); 
legend('R1','R2','BFP','RFP','GFP')

figure()
subplot(3,2,1)
plot(T_RFP,Y_RFP(:,1),'m',T_RFP,Y_RFP(:,4),'y','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOT1'); legend('R1','R4')

subplot(3,2,2)
plot(T_RFP,Y_RFP(:,2),'c',T_RFP,Y_RFP(:,3),'y','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOT2'); legend('R2','R3')

subplot(3,2,3)
plot(T_RFP,Y_RFP(:,1),'m',T_RFP,Y_RFP(:,3),'y',T_RFP,Y_RFP(:,6),'b',T_RFP,upper_threshold_RFP,'k--','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOR1: RFP=NOR(R1,R3)'); legend('R1','R3','BFP')

subplot(3,2,4)
plot(T_RFP,Y_RFP(:,2),'c',T_RFP,Y_RFP(:,4),'y',T_RFP,Y_RFP(:,7),'r',T_RFP,upper_threshold_RFP,'k--','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOR2: RFP=NOR(R2,R4)'); legend('R2','R4','RFP')

subplot(3,2,5)
plot(T_RFP,Y_RFP(:,5),'k',T_RFP,Y_RFP(:,8),'g',T_RFP,upper_threshold_RFP,'r--','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOT3: GFP=NOT(R5)'); legend('R5','GFP')

subplot(3,2,6)
plot(T_RFP,Y_RFP(:,1),'m',T_RFP,Y_RFP(:,2),'c',T_RFP,Y_RFP(:,6),'b',T_RFP,Y_RFP(:,7),'r',T_RFP,Y_RFP(:,8),'g',T_RFP,upper_threshold_RFP,'k--')
xlabel('Minutes'); ylabel('Molecules/Cell');
title('Comparator: R:R1<R2 B:R1>R2 G:R1=R2'); 
legend('R1','R2','BFP','RFP','GFP')

figure()
subplot(3,2,1)
plot(T_GFP,Y_GFP(:,1),'m',T_GFP,Y_GFP(:,4),'y','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOT1'); legend('R1','R4')

subplot(3,2,2)
plot(T_GFP,Y_GFP(:,2),'c',T_GFP,Y_GFP(:,3),'y','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOT2'); legend('R2','R3')

subplot(3,2,3)
plot(T_GFP,Y_GFP(:,1),'m',T_GFP,Y_GFP(:,3),'y',T_GFP,Y_GFP(:,6),'b',T_GFP,upper_threshold_GFP,'k--','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOR1: BFP=NOR(R1,R3)'); legend('R1','R3','BFP')

subplot(3,2,4)
plot(T_GFP,Y_GFP(:,2),'c',T_GFP,Y_GFP(:,4),'y',T_GFP,Y_GFP(:,7),'r',T_GFP,upper_threshold_GFP,'k--','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOR2: RFP=NOR(R2,R4)'); legend('R2','R4','RFP')

subplot(3,2,5)
plot(T_GFP,Y_GFP(:,5),'k',T_GFP,Y_GFP(:,8),'g',T_GFP,upper_threshold_GFP,'r--','LineWidth',2)
xlabel('Minutes'); ylabel('Molecules/Cell');
title('NOT3: GFP=NOT(R5)'); legend('R5','GFP')

subplot(3,2,6)
plot(T_GFP,Y_GFP(:,1),'m',T_GFP,Y_GFP(:,2),'c',T_GFP,Y_GFP(:,6),'b',T_GFP,Y_GFP(:,7),'r',T_GFP,Y_GFP(:,8),'g',T_GFP,upper_threshold_GFP,'k--')
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