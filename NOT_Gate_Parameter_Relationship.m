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
k_R5 = 300*k_dR5;
k_BFP = 300*k_dBFP;
k_RFP = 300*k_dRFP;
k_GFP = 300*k_dGFP;

K_R1 = 250;
K_R2 = 250;
K_R1R3 = K_R1;
K_R3R1 = 250;
K_R2R4 = K_R2;
K_R4R2 = 250;
K_R5 = 250;

n_R1 = 3;
n_R2 = 3;
n_R3 = 3;
n_R4 = 3;
n_R5 = 3;

%GATE-PARAMETER RELATIONSHIP TRANSFER FXNS
input_range = 1000;
input_step = 50;
parameter_start = 100;
parameter_range = 500;
parameter_step = 50;
gate_output = zeros(length(0:input_step:input_range),length(parameter_start:parameter_step:parameter_range));
gate_input = zeros(length(0:input_step:input_range),length(parameter_start:parameter_step:parameter_range));

lower_threshold = 100 * ones(length(gate_input));
upper_threshold = 400 * ones(length(gate_input));
lower_threshold_vector = linspace(0,300,length(gate_input)); 
upper_threshold_vector = linspace(0,300,length(gate_input));

%NOT GATES 1/2
% m = 0;
% for parameter = parameter_start:parameter_step:parameter_range
%     m = m+1;
%     n = 0;
%     for input_species = 0:input_step:input_range
%         n = n+1;
%          
%         not1=BioSystem();
%         not1.AddConstant('k_dR2', k_dR2);
%         not1.AddConstant('k_R3', k_R3);
%         not1.AddConstant('k_dR3', k_dR3);
%         not1.AddConstant('K_R2',parameter);
%         not1.AddConstant('n_R2', n_R2);
% 
%         dR2dt = not1.AddCompositor('R2', input_species);
%         dR3dt = not1.AddCompositor('R3', 0);   
% 
%         not1.AddPart(Part('R2-|R3', dR3dt,...
%             Rate('k_R3*( (K_R2^n_R2)/( (K_R2^n_R2)+(R2^n_R2) ) ) - k_dR3*R3') ));
% 
%         [T,Y] = not1.run([0,500]);
%         
%         gate_input(n,m) = input_species;
%         gate_output(n,m)=Y(end,not1.CompositorIndex('R3'));
%     end
% end
% 
% figure()
% hold on
% for i = 1:1:m
%     plot(gate_input(:,i),gate_output(:,i),'LineWidth',2)
%     xlabel('molecules of R2'); ylabel('molecules of R3'); title('NOT1 Transfer Fxn')
%     legend('150','250','350','450','550','650')
% end
% plot(lower_threshold(:,1),lower_threshold_vector','r--',upper_threshold(:,1),upper_threshold_vector','r--','LineWidth',2)

% m = 0;
% for parameter = parameter_start:parameter_step:parameter_range
%     m = m+1;
%     n = 0;
%     for input_species = 0:input_step:input_range
%         n = n+1;
%          
%         not3=BioSystem();
%         not3.AddConstant('k_dR5', k_dR5);
%         not3.AddConstant('k_GFP', k_GFP);
%         not3.AddConstant('k_dGFP', k_dGFP);
%         not3.AddConstant('K_R5',parameter);
%         not3.AddConstant('n_R5', n_R5);
% 
%         dR5dt = not3.AddCompositor('R5', input_species);
%         dGFPdt = not3.AddCompositor('GFP', 0);   
% 
%         not3.AddPart(Part('R5-|GFP', dGFPdt,...
%             Rate('k_GFP*( (K_R5^n_R5)/( (K_R5^n_R5)+(R5^n_R5) )* ) - k_dGFP*GFP') ));
% 
%         [T,Y] = not3.run([0,500]);
%         
%         gate_input(n,m) = input_species;
%         gate_output(n,m)=Y(end,not3.CompositorIndex('GFP'));
%     end
% end


%Simulating additional tandem operator sites
parameter_start = 1;
parameter_step = 1;
parameter_range = 5;

m = 0;
for parameter = parameter_start:parameter_step:parameter_range
    m = m+1;
    n = 0;
    for input_species = 0:input_step:input_range
        n = n+1;
         
        not3=BioSystem();
        not3.AddConstant('k_dR5', k_dR5);
        not3.AddConstant('k_GFP', k_GFP);
        not3.AddConstant('k_dGFP', k_dGFP);
        not3.AddConstant('K_R5',K_R5);
        not3.AddConstant('n_R5', n_R5);
        not3.AddConstant('n_tandem',parameter);

        dR5dt = not3.AddCompositor('R5', input_species);
        dGFPdt = not3.AddCompositor('GFP', 0);   

        not3.AddPart(Part('R5-|GFP', dGFPdt,...
            Rate('k_GFP*( (K_R5^n_R5)/( (K_R5^n_R5)+(R5^n_R5) ) )^(n_tandem) - k_dGFP*GFP') ));

        [T,Y] = not3.run([0,500]);
        
        gate_input(n,m) = input_species;
        gate_output(n,m)=Y(end,not3.CompositorIndex('GFP'));
    end
end

figure()
hold on
for i = 1:1:m
    plot(gate_input(:,i),gate_output(:,i),'LineWidth',2)
    xlabel('molecules of R5'); ylabel('molecules of GFP'); title('NOT3 Transfer Fxn')
    legend('1','2','3','4','5')
end
plot(lower_threshold(:,1),lower_threshold_vector','r--',upper_threshold(:,1),upper_threshold_vector','r--','LineWidth',2)

%figure(3)
%loglog(gate_input(:,1),gate_output(:,1),gate_input(:,2),gate_output(:,2)) 

% figure()
% subplot(2,1,1)
% plot(T,Y(:,not1.CompositorIndex('R1')),'r*',T, Y(:,not1.CompositorIndex('R3')),'b*')
% xlabel('minutes'); ylabel('molecules'); legend('R1','R2')
% 
% subplot(2,1,2)
% plot(input,output,'b--','LineWidth',2) 
% xlabel('[R1]'); ylabel('[R3]'); title('Not1 Transfer Fxn')
