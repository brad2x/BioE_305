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
input_range_1 = 200;
input_step_1 = 40;

input_start_2 = 500;
input_range_2 = 1500;
input_step_2 = 200;

BFP_output = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));
RFP_output = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));
GFP_output = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));
R5_output = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));
R3_output = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));
R4_output = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));

x_scaling = input_start_1:input_step_1:input_range_1;   
y_scaling = input_start_2:input_step_2:input_range_2;

m =0;
for input_species_1 = input_start_1:input_step_1:input_range_1
    m = m+1;
    n=0;
    for input_species_2 = input_start_2:input_step_2:input_range_2
        n = n+1;

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

        dR1dt = comparator.AddCompositor('R1', input_species_1);
        dR2dt = comparator.AddCompositor('R2', input_species_2);
        dR3dt = comparator.AddCompositor('R3', 0);
        dR4dt = comparator.AddCompositor('R4', 0);
        dR5dt = comparator.AddCompositor('R5', 0);
        dBFPdt = comparator.AddCompositor('BFP', 0);
        dRFPdt = comparator.AddCompositor('RFP', 0);
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
        
        [T,Y] = comparator.run([0,500]);
        BFP_output(n,m)=Y(end,comparator.CompositorIndex('BFP'));
        RFP_output(n,m)=Y(end,comparator.CompositorIndex('RFP'));
        GFP_output(n,m)=Y(end,comparator.CompositorIndex('GFP'));
        R5_output(n,m)=Y(end,comparator.CompositorIndex('R5'));
        R3_output(n,m)=Y(end,comparator.CompositorIndex('R3'));
        R4_output(n,m)=Y(end,comparator.CompositorIndex('R4'));
    end
end

% figure()
% ax1 = subplot(3,2,1);
% contourf(x_scaling,y_scaling,BFP_output)
% map1 = [0, 0, 0
%     0, 0, 0.05
%     0, 0, 0.1
%     0, 0, 0.15
%     0, 0, 0.2
%     0, 0, 0.25
%     0, 0, 0.3
%     0, 0, 0.35
%     0, 0, 0.4
%     0, 0, 0.45
%     0, 0, 0.5
%     0, 0, 0.55
%     0, 0, 0.6
%     0, 0, 0.65
%     0, 0, 0.7
%     0, 0, 0.75
%     0, 0, 0.8
%     0, 0, 0.85
%     0, 0, 0.9
%     0, 0, 0.95
%     0, 0, 1.0];
% colormap(ax1, map1)
% colorbar
% xlabel('R1 molecules/cell'); ylabel('R2 molecules/cell'); title('BFP Output (mlcl/cell)')
% 
% ax2 = subplot(3,2,2);
% contourf(x_scaling,y_scaling,RFP_output)
% map2 = [0, 0, 0
%     0.05, 0, 0
%     0.1, 0, 0
%     0.15, 0, 0
%     0.2, 0, 0
%     0.25, 0, 0
%     0.3, 0, 0
%     0.35, 0, 0
%     0.4, 0, 0
%     0.45, 0, 0
%     0.5, 0, 0
%     0.55, 0, 0
%     0.6, 0, 0
%     0.65, 0, 0
%     0.7, 0, 0
%     0.75, 0, 0
%     0.8, 0, 0
%     0.85, 0, 0
%     0.9, 0, 0
%     0.95, 0, 0
%     1.0, 0, 0];
% colormap(ax2, map2)
% colorbar
% xlabel('R1 molecules/cell'); ylabel('R2 molecules/cell'); title('RFP Output (mlcl/cell)')
% 
% ax3 = subplot(3,2,4);
% contourf(x_scaling,y_scaling,GFP_output)
% map3 = [0, 0, 0
%     0, 0.05, 0
%     0, 0.1, 0
%     0, 0.15, 0
%     0, 0.2, 0
%     0, 0.25, 0
%     0, 0.3, 0
%     0, 0.35, 0
%     0, 0.4, 0
%     0, 0.45, 0
%     0, 0.5, 0
%     0, 0.55, 0
%     0, 0.6, 0
%     0, 0.65, 0
%     0, 0.7, 0
%     0, 0.75, 0
%     0, 0.8, 0
%     0, 0.85, 0
%     0, 0.9, 0
%     0, 0.95, 0
%     0, 1.0, 0];
% colormap(ax3, map3)
% colorbar
% xlabel('R1 molecules/cell'); ylabel('R2 molecules/cell'); title('GFP Output (mlcl/cell)')
% 
% ax4 = subplot(3,2,3);
% contourf(x_scaling,y_scaling,R5_output)
% map4 = [0, 0, 0
%     0.05, 0, 0.05
%     0.1, 0, 0.1
%     0.15, 0, 0.15
%     0.2, 0, 0.2
%     0.25, 0, 0.25
%     0.3, 0, 0.3
%     0.35, 0, 0.35
%     0.4, 0, 0.4
%     0.45, 0, 0.45
%     0.5, 0, 0.5
%     0.55, 0, 0.55
%     0.6, 0, 0.6
%     0.65, 0, 0.65
%     0.7, 0, 0.7
%     0.75, 0, 0.75
%     0.8, 0, 0.8
%     0.85, 0, 0.85
%     0.9, 0, 0.9
%     0.95, 0, 0.95
%     1.0, 0, 1.0];
% colormap(ax4, map4)
% colorbar
% xlabel('R1 molecules/cell'); ylabel('R2 molecules/cell'); title('R5 (mlcl/cell)')
% 
% ax4 = subplot(3,2,5);
% contourf(x_scaling,y_scaling,R3_output)
% colormap(ax4, gray)
% colorbar
% xlabel('R1 molecules/cell'); ylabel('R2 molecules/cell'); title('R3 (mlcl/cell)')
% 
% ax4 = subplot(3,2,6);
% contourf(x_scaling,y_scaling,R4_output)
% colormap(ax4, gray)
% colorbar
% xlabel('R1 molecules/cell'); ylabel('R2 molecules/cell'); title('R4 (mlcl/cell)')

GFP_threshold_high = 290;
GFP_threshold_low = 290;
RFP_threshold = 290;
BFP_threshold = 290;

for i = 1:1:length(BFP_output(:,1))
    for k = 1:1:length(BFP_output(1,:))
        if BFP_output(i,k) >= BFP_threshold
            %BFP_output(i,k) = 1;
        elseif BFP_output(i,k) < BFP_threshold
            BFP_output(i,k) = 0;
        end
        
        if RFP_output(i,k) >= RFP_threshold
            %RFP_output(i,k) = 1;
        elseif RFP_output(i,k) < RFP_threshold
            RFP_output(i,k) = 0;
        end
 
        if GFP_output_high(i,k) >= GFP_threshold_high
            %GFP_output_high(i,k) = 1;
        elseif GFP_output_high(i,k) < GFP_threshold_high
            GFP_output_high(i,k) = 0;
        end
        
        if GFP_output_low(i,k) >= GFP_threshold_low
            %GFP_output_low(i,k) = 1;
        elseif GFP_output_low(i,k) < GFP_threshold_low
            GFP_output_low(i,k) = 0;
        end
    end
end

figure()
contourf(x_scaling,y_scaling,BFP_output)
map5 = [1, 1, 1
    0, 0, 1];
colormap(map5)
colorbar
xlabel('R1 molecules/cell'); ylabel('R2 molecules/cell'); title('BFP Output (mlcl/cell)')

figure()
contourf(x_scaling,y_scaling,RFP_output)
map6 = [1, 1, 1
    1, 0, 0];
colormap(map6)
colorbar
xlabel('R1 molecules/cell'); ylabel('R2 molecules/cell'); title('RFP Output (mlcl/cell)')

figure()
contourf(x_scaling,y_scaling,GFP_output)
map7 = [1, 1, 1
    0, 1, 0];
colormap(map7)
colorbar
xlabel('R1 molecules/cell'); ylabel('R2 molecules/cell'); title('GFP Output (mlcl/cell)')