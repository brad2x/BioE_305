urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
         'part_composition_setup.m');
rehash;
part_composition_setup('v5');

k_dR1 = .0277;
k_dR2 = .0277;
k_dR3 = .0277;
k_dR4 = .0277;
k_dR5 = .0277;
k_dBFP = .0277;
k_dRFP = .0277;
k_dGFP = .0277;

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

input_start_1 = 0;
input_range_1 = 600;
input_step_1 = 25;

input_2_low = 0;
input_2_high = 1000;

parameter_start = 100;
parameter_step = 50;
parameter_range = 500;

gate_input = zeros(length(input_start_1:input_step_1:input_range_1),length(parameter_start:parameter_step:parameter_range));

nor1_gate_output_low = zeros(length(input_start_1:input_step_1:input_range_1),length(parameter_start:parameter_step:parameter_range));
nor1_gate_output_high = zeros(length(input_start_1:input_step_1:input_range_1),length(parameter_start:parameter_step:parameter_range));
nor2_gate_output_low = zeros(length(input_start_1:input_step_1:input_range_1),length(parameter_start:parameter_step:parameter_range));
nor2_gate_output_high = zeros(length(input_start_1:input_step_1:input_range_1),length(parameter_start:parameter_step:parameter_range));

lower_threshold = 50 * ones(length(gate_input));
upper_threshold = 450 * ones(length(gate_input));
lower_threshold_vector = linspace(0,500,length(gate_input)); 
upper_threshold_vector = linspace(0,500,length(gate_input));

m = 0;
for parameter = parameter_start:parameter_step:parameter_range
    m = m+1;
    n=0;
    for input_species_1 = input_start_1:input_step_1:input_range_1
        n = n+1;

        nor1l=BioSystem();
        nor1l.AddConstant('k_dR1', k_dR1);
        nor1l.AddConstant('k_dR3', k_dR3);
        nor1l.AddConstant('k_R5', k_R5);
        nor1l.AddConstant('k_dR5', k_dR5);
        nor1l.AddConstant('K_R1R3',K_R1);
        nor1l.AddConstant('n_R1', n_R1);
        nor1l.AddConstant('K_R3R1',parameter);
        nor1l.AddConstant('n_R3', n_R3);

        dR1dt = nor1l.AddCompositor('R1', input_2_low);
        dR3dt = nor1l.AddCompositor('R3', input_species_1);
        dR5dt = nor1l.AddCompositor('R5', 0);   

        nor1l.AddPart(Part('R1,R3-|R5', dR5dt,...
            Rate('k_R5*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) ) - k_dR5*R5') ));

        [T,Y] = nor1l.run([0,200]);
        gate_input(n,m) = input_species_1;
        nor1_gate_output_low(n,m)=Y(end,nor1l.CompositorIndex('R5'));
    end
end

m = 0;
for parameter = parameter_start:parameter_step:parameter_range
    m = m+1;
    n=0;
    for input_species_1 = input_start_1:input_step_1:input_range_1
        n = n+1;

        nor1h=BioSystem();
        nor1h.AddConstant('k_dR1', k_dR1);
        nor1h.AddConstant('k_dR3', k_dR3);
        nor1h.AddConstant('k_R5', k_R5);
        nor1h.AddConstant('k_dR5', k_dR5);
        nor1h.AddConstant('K_R1R3',K_R1);
        nor1h.AddConstant('n_R1', n_R1);
        nor1h.AddConstant('K_R3R1',parameter);
        nor1h.AddConstant('n_R3', n_R3);

        dR1dt = nor1h.AddCompositor('R1', input_2_high);
        dR3dt = nor1h.AddCompositor('R3', input_species_1);
        dR5dt = nor1h.AddCompositor('R5', 0);   

        nor1h.AddPart(Part('R1,R3-|R5', dR5dt,...
            Rate('k_R5*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) ) - k_dR5*R5') ));

        [T,Y] = nor1h.run([0,200]);
        nor1_gate_output_high(n,m)=Y(end,nor1h.CompositorIndex('R5'));
    end
end

figure()
%subplot(2,1,1)
hold on
for i = 1:1:m
    plot(gate_input(:,i),nor1_gate_output_low(:,i))
    xlabel('molecules of R3'); ylabel('molecules of R5');
    title('NOR1 at [R1]=100'); legend('100','150','200','250','300','350','400','450','500');
end
plot(lower_threshold(:,1),lower_threshold_vector','r--',upper_threshold(:,1),upper_threshold_vector','r--','LineWidth',2)

figure()
%subplot(2,1,2)
hold on
for i = 1:1:m
    plot(gate_input(:,i),nor1_gate_output_high(:,i))
    xlabel('molecules of R3'); ylabel('molecules of R5');
    title('NOR1 at [R1]=1000'); legend('100','150','200','250','300','350','400','450','500');
end
plot(lower_threshold(:,1),.05*lower_threshold_vector','r--',upper_threshold(:,1),.05*upper_threshold_vector','r--','LineWidth',2)
