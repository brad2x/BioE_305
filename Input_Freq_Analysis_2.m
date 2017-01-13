urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
         'part_composition_setup.m');
rehash;
part_composition_setup('v5');

k_dR1 = log(2)/25;
k_dR2 = log(2)/25;
k_dR3 = log(2)/5;
k_dR4 = log(2)/5;
k_dR5 = log(2)/5;
k_dBFP = log(2)/5;
k_dRFP = log(2)/5;
k_dGFP = log(2)/5;

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

parameter_start = 40;
parameter_range = 100;
parameter_step = 4;

parameter_input = zeros(length(parameter_start:parameter_step:parameter_range),1);
pass_threshold = zeros(length(parameter_start:parameter_step:parameter_range),1);

n = 0;
for parameter = parameter_start:parameter_step:parameter_range
    n = n+1;
    
    A2 = 900/2;
    w2 = (2*pi)/(parameter);
    
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

    %comparator.AddConstant('A1', A1);
    %comparator.AddConstant('w1', w1);
    comparator.AddConstant('A2', A2);
    comparator.AddConstant('w2', w2);

    dR1dt = comparator.AddCompositor('R1', 1000);
    dR2dt = comparator.AddCompositor('R2', 100);
    dR3dt = comparator.AddCompositor('R3', 0);
    dR4dt = comparator.AddCompositor('R4', 0);
    dR5dt = comparator.AddCompositor('R5', 0);
    dBFPdt = comparator.AddCompositor('BFP', 0);
    dRFPdt = comparator.AddCompositor('RFP', 0);
    dGFPdt = comparator.AddCompositor('GFP', 0);

    dtimedt = comparator.AddCompositor('T1',0);

    comparator.AddPart(Part('time', dtimedt, Rate('1')));

    %comparator.AddPart(Part('R1', dR1dt,Rate('-A1*w1*sin(w1*T1)')));

    comparator.AddPart(Part('R2', dR2dt,Rate('-A2*w2*sin(w2*T1 + pi)')));

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

    [T,Y] = comparator.run([0,1.5*parameter_range]);

    threshold = 275*ones(length(T),1);
    
    parameter_input(n) = parameter;
    
%     if max(Y(:,comparator.CompositorIndex('GFP'))) >= threshold(1)
%         pass_threshold(n) = 1;
%     else
%         pass_threshold(n) = 0;
%     end
    
    min_above_threshold = find(Y(:,comparator.CompositorIndex('GFP')) >= threshold(1),1);
    if( isempty(min_above_threshold) )
        pass_threshold(n) = parameter_range;
    elseif( (min_above_threshold - .5*parameter) > (1.5*parameter))
        pass_threshold(n) = parameter_range;
    else
        pass_threshold(n) = T(min_above_threshold) - .5*parameter;
    end
    
end 

figure()
plot(parameter_input, pass_threshold,'*')
xlabel('Oscillation Period'); ylabel('GFP Delay');
title('Input Freq vs Threshold');

% figure()
% plot(T,Y(:,1),T,Y(:,2),T,Y(:,comparator.CompositorIndex('GFP')),'g',T,threshold,'k--')
% legend('R1','R2','GFP')