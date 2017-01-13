urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
         'part_composition_setup.m');
rehash;
part_composition_setup('v5');

k_R1 = 13.85;
k_R2 = 13.85;
k_R3 = 13.85;
k_R4 = 13.85;
k_R5 = 13.85;
k_BFP = 13.85;
k_RFP = 13.85;
k_GFP = 13.85;

k_dR1 = .0277;
k_dR2 = .0277;
k_dR3 = .0277;
k_dR4 = .0277;
k_dR5 = .0277;
k_dBFP = .0277;
k_dRFP = .0277;
k_dGFP = .0277;

K_R1 = 250;
K_R2 = 250;
K_R3 = 250;
K_R4 = 250;
K_R1R3 = 250;
K_R3R1 = 250;
K_R2R4 = 250;
K_R4R2 = 250;
K_R5 = 250;

n_R1 = 3;
n_R2 = 3;
n_R3 = 3;
n_R4 = 3;
n_R5 = 3;

input_start_1 = 0;
input_range_1 = 1000;
input_step_1 = 100;

input_start_2 = 0;
input_range_2 = 1000;
input_step_2 = 100;


nor1_gate_output = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));
x_scaling = input_start_1:input_step_1:input_range_1;   
y_scaling = input_start_2:input_step_2:input_range_2;

lower_threshold = 100 * ones(length(gate_input));
upper_threshold = 400 * ones(length(gate_input));
lower_threshold_vector = linspace(0,150,length(gate_input)); 
upper_threshold_vector = linspace(0,150,length(gate_input));

m =0;
for input_species_1 = input_start_1:input_step_1:input_range_1
    m = m+1;

    n=0;
    for input_species_2 = input_start_2:input_step_2:input_range_2
        n = n+1;

        nor1=BioSystem();
        nor1.AddConstant('k_dR1', k_dR1);
        nor1.AddConstant('k_dR3', k_dR3);
        nor1.AddConstant('k_R5', k_R5);
        nor1.AddConstant('k_dR5', k_dR5);
        nor1.AddConstant('K_R1R3',K_R1);
        nor1.AddConstant('n_R1', n_R1);
        nor1.AddConstant('K_R3R1',K_R3);
        nor1.AddConstant('n_R3', n_R3);
        nor1.AddConstant('k_bR5', .277);    %Basal production set to produce 10 mlcl at s-s

        dR1dt = nor1.AddCompositor('R1', input_species_1);
        dR3dt = nor1.AddCompositor('R3', input_species_2);
        dR5dt = nor1.AddCompositor('R5', 0);   

        nor1.AddPart(Part('R1,R3-|R5', dR5dt,...
            Rate('k_R5*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) ) - k_dR5*R5 + k_bR5') ));

        [T,Y] = nor1.run([0,100]);
        nor1_gate_output(n,m)=Y(end,nor1.CompositorIndex('R5'));
    end
end

figure()
contourf(x_scaling,y_scaling,nor1_gate_output)
colorbar
xlabel('[R1]'); ylabel('[R3]'); title('Nor1 Transfer Fxn')