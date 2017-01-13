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

% ALL GATE TRANSFER FUNCTIONS
input_range = 1200;
input_step = 100;

num_gates = 3;

gate_output = zeros(length(0:input_step:input_range), num_gates);
gate_input = zeros(length(0:input_step:input_range), num_gates);

lower_threshold = 100 * ones(length(gate_input));
upper_threshold = 1000 * ones(length(gate_input));
lower_threshold_vector = linspace(0,300,length(gate_input)); 
upper_threshold_vector = linspace(0,300,length(gate_input));

% g=2;
% n=0;
% for input_species = 0:input_step:input_range
%     n = n+1;
% 
%     not2=BioSystem();
%     not2.AddConstant('k_dR1', k_dR1);
%     not2.AddConstant('k_R4', k_R4);
%     not2.AddConstant('k_dR4', k_dR4);
%     not2.AddConstant('K_R1',K_R1);
%     not2.AddConstant('n_R1', n_R1);
%     not2.AddConstant('k_bR4', 0);    %Basal production set to produce 10 mlcl at s-s
% 
%     dR1dt = not2.AddCompositor('R1', input_species);
%     dR4dt = not2.AddCompositor('R4', 0);   
% 
%     not2.AddPart(Part('R1-|R4', dR4dt,...
%         Rate('k_R4*( (K_R1^n_R1)/( (K_R1^n_R1)+(R1^n_R1) ) ) - k_dR4*R4 + k_bR4') ));
% 
%     [T,Y] = not2.run([0,250]);
%     gate_input(n,g) = input_species;
%     gate_output(n,g)=Y(end,not2.CompositorIndex('R4'));
% end
% 
% g=1;
% n=0;
% for input_species = 0:input_step:input_range
%     n = n+1;
% 
%     not1=BioSystem();
%     not1.AddConstant('k_dR2', k_dR2);
%     not1.AddConstant('k_R3', k_R3);
%     not1.AddConstant('k_dR3', k_dR3);
%     not1.AddConstant('K_R2',K_R2);
%     not1.AddConstant('n_R2', n_R2);
%     not1.AddConstant('k_bR3', 0);    %Basal production set to produce 10 mlcl at s-s
% 
%     dR2dt = not1.AddCompositor('R2', input_species);
%     dR3dt = not1.AddCompositor('R3', 0);   
% 
%     not1.AddPart(Part('R2-|R3', dR3dt,...
%         Rate('k_R3*( (K_R2^n_R2)/( (K_R2^n_R2)+(R2^n_R2) ) ) - k_dR3*R3 + k_bR3') ));
% 
%     [T,Y] = not1.run([0,250]);
%     gate_input(n,g) = input_species;
%     gate_output(n,g)=Y(end,not1.CompositorIndex('R3'));
% end

input_range = 1200;
input_step = 10;

g=3;
n=0;
for input_species = 0:input_step:input_range
    n = n+1;

    not3=BioSystem();
    not3.AddConstant('k_dR5', k_dR5);
    not3.AddConstant('k_GFP', k_GFP);
    not3.AddConstant('k_dGFP', k_dGFP);
    not3.AddConstant('K_R5',K_R5);
    not3.AddConstant('n_R5', n_R5);
    not3.AddConstant('k_bGFP', 0);    %Basal production set to produce 10 mlcl at s-s

    dR5dt = not3.AddCompositor('R5', input_species);
    dGFPdt = not3.AddCompositor('GFP', 0);   

    not3.AddPart(Part('R5-|GFP', dGFPdt,...
        Rate('k_GFP*(((K_R5^n_R5)/( (K_R5^n_R5)+(R5^n_R5)))^1) - k_dGFP*GFP + k_bGFP') ));

    [T,Y] = not3.run([0,250]);
    gate_input(n,g) = input_species;
    gate_output(n,g)=Y(end,not3.CompositorIndex('GFP'));
end

% input_start_1 = 0;
% input_range_1 = 1200;
% input_step_1 = 100;
% 
% input_start_2 = 0;
% input_range_2 = 600;
% input_step_2 = 50;
% 
% nor1_gate_output = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));
% x_scaling = input_start_1:input_step_1:input_range_1;   
% y_scaling = input_start_2:input_step_2:input_range_2;

% m =0;
% for input_species_1 = input_start_1:input_step_1:input_range_1
%     m = m+1;
% 
%     n=0;
%     for input_species_2 = input_start_2:input_step_2:input_range_2
%         n = n+1;
% 
%         nor1=BioSystem();
%         nor1.AddConstant('k_dR1', k_dR1);
%         nor1.AddConstant('k_dR3', k_dR3);
%         nor1.AddConstant('k_R5', k_R5);
%         nor1.AddConstant('k_dR5', k_dR5);
%         nor1.AddConstant('K_R1R3',K_R1R3);
%         nor1.AddConstant('n_R1', n_R1);
%         nor1.AddConstant('K_R3R1',K_R3R1);
%         nor1.AddConstant('n_R3', n_R3);
%         nor1.AddConstant('k_bR5', 0);    %Basal production set to produce 10 mlcl at s-s
% 
%         dR1dt = nor1.AddCompositor('R1', input_species_1);
%         dR3dt = nor1.AddCompositor('R3', input_species_2);
%         dR5dt = nor1.AddCompositor('R5', 0);   
% 
%         nor1.AddPart(Part('R1,R3-|R5', dR5dt,...
%             Rate('k_R5*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) ) - k_dR5*R5 + k_bR5') ));
% 
%         [T,Y] = nor1.run([0,100]);
%         nor1_gate_output(n,m)=Y(end,nor1.CompositorIndex('R5'));
%     end
% end

% nor2_gate_output = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));
% m =0;
% for input_species_1 = input_start_1:input_step_1:input_range_1
%     m = m+1;
%     n=0;
%     for input_species_2 = input_start_2:input_step_2:input_range_2
%         n = n+1;
% 
%         nor2=BioSystem();
%         nor2.AddConstant('k_dR2', k_dR2);
%         nor2.AddConstant('k_dR4', k_dR4);
%         nor2.AddConstant('k_R5', k_R5);
%         nor2.AddConstant('k_dR5', k_dR5);
%         nor2.AddConstant('K_R2R4',K_R2R4);
%         nor2.AddConstant('n_R2', n_R2);
%         nor2.AddConstant('K_R4R2',K_R4R2);
%         nor2.AddConstant('n_R4', n_R4);
%         nor2.AddConstant('k_bR5', 0);    %Basal production set to produce 10 mlcl at s-s
% 
%         dR2dt = nor2.AddCompositor('R2', input_species_1);
%         dR4dt = nor2.AddCompositor('R4', input_species_2);
%         dR5dt = nor2.AddCompositor('R5', 0);   
% 
%         nor2.AddPart(Part('R2,R4-|R5', dR5dt,...
%             Rate('k_R5*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dR5*R5 + k_bR5') ));
% 
%         [T,Y] = nor2.run([0,100]);
%         nor2_gate_output(n,m)=Y(end,nor2.CompositorIndex('R5'));
%     end
% end

% figure()
% subplot(2,3,1)
% plot(gate_input(:,1),gate_output(:,1),'b',...
%     lower_threshold,lower_threshold_vector,'r--',...
%     upper_threshold, upper_threshold_vector, 'r--','LineWidth',2) 
% xlabel('molecules of R1'); ylabel('molecules of R4'); title('NOT1 Transfer Fxn')
% 
% subplot(2,3,4)
% plot(gate_input(:,2),gate_output(:,2),'b',...
%     lower_threshold,lower_threshold_vector,'r--',...
%     upper_threshold, upper_threshold_vector, 'r--','LineWidth',2) 
% xlabel('molecules of R2'); ylabel('molecules of R3'); title('NOT2 Transfer Fxn')
% 
% subplot(2,3,3)
% plot(gate_input(:,3),gate_output(:,3),'b',...
%     'LineWidth',2) 
% xlabel('molecules of R5'); ylabel('molecules of GFP'); title('NOT3 Transfer Fxn')
% 
% subplot(2,3,2)
% contourf(x_scaling,y_scaling,nor1_gate_output)
% colorbar
% xlabel('molecules of R1'); ylabel('molecules of R3'); title('NOR1 Transfer Fxn')
% legend('[R5]')
% 
% % figure()
% subplot(2,3,5)
% contourf(x_scaling,y_scaling,nor2_gate_output)
% colorbar
% xlabel('molecules of R2'); ylabel('molecules of R4'); title('NOR2 Transfer Fxn')
% legend('[R5]')

figure()
hold on
plot(gate_input(:,3),gate_output(:,3),'b',...
    'LineWidth',2) 
plot(100*ones(length(gate_input(:,3))),linspace(0,300,length(gate_input(:,3))),'r--','LineWidth',2)
plot(400*ones(length(gate_input(:,3))),linspace(0,300,length(gate_input(:,3))),'r--','LineWidth',2)
xlabel('molecules of R5'); ylabel('molecules of GFP'); title('NOT3 Transfer Fxn')

% figure()
% contourf(x_scaling,y_scaling,nor1_gate_output)
% colorbar
% xlabel('molecules of R1'); ylabel('molecules of R3'); title('NOR1 Transfer Fxn')
% legend('[R5]')