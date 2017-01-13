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

parameter_start = 101;
parameter_range = 102;
parameter_step = .01;

parameter_input = zeros(length(parameter_start:parameter_step:parameter_range),1);
threshold_time = zeros(length(parameter_start:parameter_step:parameter_range),1);

T_store = zeros(95,length(parameter_start:parameter_step:parameter_range));
GFP_output = zeros(95,length(parameter_start:parameter_step:parameter_range));

figure()
hold on

m = 0;
for parameter = parameter_start:parameter_step:parameter_range
    m = m+1;

    comparator=BioSystem();
    comparator.AddConstant('k_R3', k_R3);
    comparator.AddConstant('k_R4', k_R4);
    comparator.AddConstant('k_R5', k_R5);
    comparator.AddConstant('k_BFP', k_BFP);
    comparator.AddConstant('k_RFP', k_RFP);
    comparator.AddConstant('k_GFP', k_GFP);

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

    dR1dt = comparator.AddCompositor('R1', 1000);
    dR2dt = comparator.AddCompositor('R2', 0);
    dR3dt = comparator.AddCompositor('R3', 500);
    dR4dt = comparator.AddCompositor('R4', 0);
    dR5dt = comparator.AddCompositor('R5', 500);
    dBFPdt = comparator.AddCompositor('BFP', 0);
    dRFPdt = comparator.AddCompositor('RFP', 300);
    dGFPdt = comparator.AddCompositor('GFP', 0);   

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

    [T,Y] = comparator.run_pulses([...
    Pulse(0,'R2',0), ...
    Pulse(100,'R2',1000), ...
    Pulse(100 + parameter,'R2',0), ...
    Pulse(400,'R2',0) ...
    %Pulse(100 + 3*parameter,'R2',0), ...
    %Pulse(100 + 4*parameter,'R2',1000), ...
    %Pulse(100 + 5*parameter,'R2',0), ...
    ]);
    
    upper_threshold = 275*ones(length(T));
    parameter_input(m) = parameter;
    
    
    for j = 1:1:length(T)
        T_store(j,m) = T(j);
        GFP_output(j,m) = Y(j,8);
    end
    
    if( max(GFP_output(:,i)) >= 275)
        color = 'g';
    else
        color = 'r';
    end
    
    plot(T,Y(:,8),color,T,upper_threshold,'k--',T,Y(:,1),'m',T,Y(:,2),'c')
    xlabel('minutes'); ylabel('number of molecules of GFP');
    legend('GFP','upper threshold','R1','R2');

    
    min_above_threshold = find(Y(:,comparator.CompositorIndex('GFP')) >= upper_threshold,1);
    if( isempty(min_above_threshold) )
        threshold_time(m) = 0;
    elseif( (min_above_threshold - 100) > (1.5*parameter))
        threshold_time(m) = 0;
    else
        threshold_time(m) = T(min_above_threshold) - 100;
    end
end

hold off

figure()
plot(parameter_input, threshold_time,'b*','LineWidth',2)
xlabel('Pulse Period'); ylabel('GFP Delay'); 
title('Frequency Response');

figure()
for i = 1:1:m
    subplot( ceil(sqrt(m)), ceil(sqrt(m)), i)
    
    if( max(GFP_output(:,i)) >= 275)
        color = 'g';
    else
        color = 'r';
    end
    
    plot(T_store(:,i),GFP_output(:,i),color,T_store(:,i),275*ones(length(T_store(:,i))),'k--')
    pulse_duration = int2str( parameter_input(i) );
    title(parameter_input(i))
    %xlabel('minutes'); ylabel('number of molecules of GFP');
    %legend('GFP','upper threshold');
    %axis([0 400 250 300]);
end