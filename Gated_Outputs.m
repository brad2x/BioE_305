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

input_start_2 = 900;
input_range_2 = 1100;
input_step_2 = 40;

BFP_output = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));

x_scaling_BFP = input_start_1:input_step_1:input_range_1;   
y_scaling_BFP = input_start_2:input_step_2:input_range_2;

m =0;
for input_species_1 = input_start_1:input_step_1:input_range_1
    m = m+1;
    n=0;
    for input_species_2 = input_start_2:input_step_2:input_range_2
        n = n+1;

        BFPcomp=BioSystem();
        
        BFPcomp.AddConstant('k_R3', k_R3);
        BFPcomp.AddConstant('k_R4', k_R4);
        BFPcomp.AddConstant('k_R5', k_R5);
        BFPcomp.AddConstant('k_BFP', k_BFP);
        BFPcomp.AddConstant('k_RFP', k_RFP);
        BFPcomp.AddConstant('k_GFP', k_GFP);

        BFPcomp.AddConstant('k_dR1', k_dR1);
        BFPcomp.AddConstant('k_dR2', k_dR2);
        BFPcomp.AddConstant('k_dR3', k_dR3);
        BFPcomp.AddConstant('k_dR4', k_dR4);
        BFPcomp.AddConstant('k_dR5', k_dR5);
        BFPcomp.AddConstant('k_dBFP', k_dBFP);
        BFPcomp.AddConstant('k_dRFP', k_dRFP);
        BFPcomp.AddConstant('k_dGFP', k_dGFP);

        BFPcomp.AddConstant('K_R1',K_R1);
        BFPcomp.AddConstant('K_R2',K_R2);
        BFPcomp.AddConstant('K_R1R3',K_R1R3);
        BFPcomp.AddConstant('K_R3R1',K_R3R1);
        BFPcomp.AddConstant('K_R2R4',K_R2R4);
        BFPcomp.AddConstant('K_R4R2',K_R4R2);
        BFPcomp.AddConstant('K_R5',K_R5);

        BFPcomp.AddConstant('n_R1', n_R1);
        BFPcomp.AddConstant('n_R2', n_R2);
        BFPcomp.AddConstant('n_R3', n_R3);
        BFPcomp.AddConstant('n_R4', n_R4);
        BFPcomp.AddConstant('n_R5', n_R5);

        dR1dt = BFPcomp.AddCompositor('R1', input_species_1);
        dR2dt = BFPcomp.AddCompositor('R2', input_species_2);
        dR3dt = BFPcomp.AddCompositor('R3', 0);
        dR4dt = BFPcomp.AddCompositor('R4', 0);
        dR5dt = BFPcomp.AddCompositor('R5', 0);
        dBFPdt = BFPcomp.AddCompositor('BFP', 0);
        dRFPdt = BFPcomp.AddCompositor('RFP', 0);
        dGFPdt = BFPcomp.AddCompositor('GFP', 0);

        BFPcomp.AddPart(Part('NOT1', dR4dt,...
            Rate('k_R4*( ((K_R1)^(n_R1))/( ((K_R1)^(n_R1))+((R1)^(n_R1)) ) ) - k_dR4*R4')));

        BFPcomp.AddPart(Part('NOT2', dR3dt,...
            Rate('k_R3*( ((K_R2)^(n_R2))/( ((K_R2)^(n_R2))+((R2)^(n_R2)) ) ) - k_dR3*R3')));

        BFPcomp.AddPart(Part('NOR1', [dR5dt, dBFPdt],...
            [Rate('k_R5*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) )'), ...
            Rate('k_BFP*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) ) - k_dBFP*BFP')]));

        BFPcomp.AddPart(Part('NOR2', [dR5dt, dRFPdt],...
            [Rate('k_R5*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dR5*R5'), ...
            Rate('k_RFP*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dRFP*RFP')]));

        BFPcomp.AddPart(Part('NOT3', dGFPdt,...
            Rate('k_GFP*(((K_R5^n_R5)/((K_R5^n_R5)+(R5^n_R5)))^2) - k_dGFP*GFP')));
        
        [T,Y] = BFPcomp.run([0,500]);
        BFP_output(n,m)=Y(end,BFPcomp.CompositorIndex('BFP'));
    end
end

input_start_1 = 900;
input_range_1 = 1100;
input_step_1 = 40;

input_start_2 = 0;
input_range_2 = 200;
input_step_2 = 40;

RFP_output = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));

x_scaling_RFP = input_start_1:input_step_1:input_range_1;   
y_scaling_RFP = input_start_2:input_step_2:input_range_2;

m =0;
for input_species_1 = input_start_1:input_step_1:input_range_1
    m = m+1;
    n=0;
    for input_species_2 = input_start_2:input_step_2:input_range_2
        n = n+1;

        RFPcomp=BioSystem();
        
        RFPcomp.AddConstant('k_R3', k_R3);
        RFPcomp.AddConstant('k_R4', k_R4);
        RFPcomp.AddConstant('k_R5', k_R5);
        RFPcomp.AddConstant('k_BFP', k_BFP);
        RFPcomp.AddConstant('k_RFP', k_RFP);
        RFPcomp.AddConstant('k_GFP', k_GFP);

        RFPcomp.AddConstant('k_dR1', k_dR1);
        RFPcomp.AddConstant('k_dR2', k_dR2);
        RFPcomp.AddConstant('k_dR3', k_dR3);
        RFPcomp.AddConstant('k_dR4', k_dR4);
        RFPcomp.AddConstant('k_dR5', k_dR5);
        RFPcomp.AddConstant('k_dBFP', k_dBFP);
        RFPcomp.AddConstant('k_dRFP', k_dRFP);
        RFPcomp.AddConstant('k_dGFP', k_dGFP);

        RFPcomp.AddConstant('K_R1',K_R1);
        RFPcomp.AddConstant('K_R2',K_R2);
        RFPcomp.AddConstant('K_R1R3',K_R1R3);
        RFPcomp.AddConstant('K_R3R1',K_R3R1);
        RFPcomp.AddConstant('K_R2R4',K_R2R4);
        RFPcomp.AddConstant('K_R4R2',K_R4R2);
        RFPcomp.AddConstant('K_R5',K_R5);

        RFPcomp.AddConstant('n_R1', n_R1);
        RFPcomp.AddConstant('n_R2', n_R2);
        RFPcomp.AddConstant('n_R3', n_R3);
        RFPcomp.AddConstant('n_R4', n_R4);
        RFPcomp.AddConstant('n_R5', n_R5);

        dR1dt = RFPcomp.AddCompositor('R1', input_species_1);
        dR2dt = RFPcomp.AddCompositor('R2', input_species_2);
        dR3dt = RFPcomp.AddCompositor('R3', 0);
        dR4dt = RFPcomp.AddCompositor('R4', 0);
        dR5dt = RFPcomp.AddCompositor('R5', 0);
        dBFPdt = RFPcomp.AddCompositor('BFP', 0);
        dRFPdt = RFPcomp.AddCompositor('RFP', 0);
        dGFPdt = RFPcomp.AddCompositor('GFP', 0);

        RFPcomp.AddPart(Part('NOT1', dR4dt,...
            Rate('k_R4*( ((K_R1)^(n_R1))/( ((K_R1)^(n_R1))+((R1)^(n_R1)) ) ) - k_dR4*R4')));

        RFPcomp.AddPart(Part('NOT2', dR3dt,...
            Rate('k_R3*( ((K_R2)^(n_R2))/( ((K_R2)^(n_R2))+((R2)^(n_R2)) ) ) - k_dR3*R3')));

        RFPcomp.AddPart(Part('NOR1', [dR5dt, dBFPdt],...
            [Rate('k_R5*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) )'), ...
            Rate('k_BFP*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) ) - k_dBFP*BFP')]));

        RFPcomp.AddPart(Part('NOR2', [dR5dt, dRFPdt],...
            [Rate('k_R5*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dR5*R5'), ...
            Rate('k_RFP*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dRFP*RFP')]));

        RFPcomp.AddPart(Part('NOT3', dGFPdt,...
            Rate('k_GFP*(((K_R5^n_R5)/((K_R5^n_R5)+(R5^n_R5)))^2) - k_dGFP*GFP')));
        
        [T,Y] = RFPcomp.run([0,500]);
        RFP_output(n,m)=Y(end,RFPcomp.CompositorIndex('RFP'));
    end
end

input_start_1 = 0;
input_range_1 = 200;
input_step_1 = 40;

input_start_2 = 0;
input_range_2 = 200;
input_step_2 = 40;

GFP_output_low = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));

x_scaling_GFPl = input_start_1:input_step_1:input_range_1;   
y_scaling_GFPl = input_start_2:input_step_2:input_range_2;

m =0;
for input_species_1 = input_start_1:input_step_1:input_range_1
    m = m+1;
    n=0;
    for input_species_2 = input_start_2:input_step_2:input_range_2
        n = n+1;

        GFPlcomp=BioSystem();
        
        GFPlcomp.AddConstant('k_R3', k_R3);
        GFPlcomp.AddConstant('k_R4', k_R4);
        GFPlcomp.AddConstant('k_R5', k_R5);
        GFPlcomp.AddConstant('k_BFP', k_BFP);
        GFPlcomp.AddConstant('k_RFP', k_RFP);
        GFPlcomp.AddConstant('k_GFP', k_GFP);

        GFPlcomp.AddConstant('k_dR1', k_dR1);
        GFPlcomp.AddConstant('k_dR2', k_dR2);
        GFPlcomp.AddConstant('k_dR3', k_dR3);
        GFPlcomp.AddConstant('k_dR4', k_dR4);
        GFPlcomp.AddConstant('k_dR5', k_dR5);
        GFPlcomp.AddConstant('k_dBFP', k_dBFP);
        GFPlcomp.AddConstant('k_dRFP', k_dRFP);
        GFPlcomp.AddConstant('k_dGFP', k_dGFP);

        GFPlcomp.AddConstant('K_R1',K_R1);
        GFPlcomp.AddConstant('K_R2',K_R2);
        GFPlcomp.AddConstant('K_R1R3',K_R1R3);
        GFPlcomp.AddConstant('K_R3R1',K_R3R1);
        GFPlcomp.AddConstant('K_R2R4',K_R2R4);
        GFPlcomp.AddConstant('K_R4R2',K_R4R2);
        GFPlcomp.AddConstant('K_R5',K_R5);

        GFPlcomp.AddConstant('n_R1', n_R1);
        GFPlcomp.AddConstant('n_R2', n_R2);
        GFPlcomp.AddConstant('n_R3', n_R3);
        GFPlcomp.AddConstant('n_R4', n_R4);
        GFPlcomp.AddConstant('n_R5', n_R5);

        dR1dt = GFPlcomp.AddCompositor('R1', input_species_1);
        dR2dt = GFPlcomp.AddCompositor('R2', input_species_2);
        dR3dt = GFPlcomp.AddCompositor('R3', 0);
        dR4dt = GFPlcomp.AddCompositor('R4', 0);
        dR5dt = GFPlcomp.AddCompositor('R5', 0);
        dBFPdt = GFPlcomp.AddCompositor('BFP', 0);
        dRFPdt = GFPlcomp.AddCompositor('RFP', 0);
        dGFPdt = GFPlcomp.AddCompositor('GFP', 0);

        GFPlcomp.AddPart(Part('NOT1', dR4dt,...
            Rate('k_R4*( ((K_R1)^(n_R1))/( ((K_R1)^(n_R1))+((R1)^(n_R1)) ) ) - k_dR4*R4')));

        GFPlcomp.AddPart(Part('NOT2', dR3dt,...
            Rate('k_R3*( ((K_R2)^(n_R2))/( ((K_R2)^(n_R2))+((R2)^(n_R2)) ) ) - k_dR3*R3')));

        GFPlcomp.AddPart(Part('NOR1', [dR5dt, dBFPdt],...
            [Rate('k_R5*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) )'), ...
            Rate('k_BFP*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) ) - k_dBFP*BFP')]));

        GFPlcomp.AddPart(Part('NOR2', [dR5dt, dRFPdt],...
            [Rate('k_R5*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dR5*R5'), ...
            Rate('k_RFP*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dRFP*RFP')]));

        GFPlcomp.AddPart(Part('NOT3', dGFPdt,...
            Rate('k_GFP*(((K_R5^n_R5)/((K_R5^n_R5)+(R5^n_R5)))^2) - k_dGFP*GFP')));
        
        [T,Y] = GFPlcomp.run([0,500]);
        GFP_output_low(n,m)=Y(end,GFPlcomp.CompositorIndex('GFP'));
    end
end

input_start_1 = 900;
input_range_1 = 1100;
input_step_1 = 40;

input_start_2 = 900;
input_range_2 = 1100;
input_step_2 = 40;

GFP_output_high = zeros(length(input_start_1:input_step_1:input_range_1), length(input_start_2:input_step_2:input_range_2));

x_scaling_GFPh = input_start_1:input_step_1:input_range_1;   
y_scaling_GFPh = input_start_2:input_step_2:input_range_2;

m =0;
for input_species_1 = input_start_1:input_step_1:input_range_1
    m = m+1;
    n=0;
    for input_species_2 = input_start_2:input_step_2:input_range_2
        n = n+1;

        GFPhcomp=BioSystem();
        
        GFPhcomp.AddConstant('k_R3', k_R3);
        GFPhcomp.AddConstant('k_R4', k_R4);
        GFPhcomp.AddConstant('k_R5', k_R5);
        GFPhcomp.AddConstant('k_BFP', k_BFP);
        GFPhcomp.AddConstant('k_RFP', k_RFP);
        GFPhcomp.AddConstant('k_GFP', k_GFP);

        GFPhcomp.AddConstant('k_dR1', k_dR1);
        GFPhcomp.AddConstant('k_dR2', k_dR2);
        GFPhcomp.AddConstant('k_dR3', k_dR3);
        GFPhcomp.AddConstant('k_dR4', k_dR4);
        GFPhcomp.AddConstant('k_dR5', k_dR5);
        GFPhcomp.AddConstant('k_dBFP', k_dBFP);
        GFPhcomp.AddConstant('k_dRFP', k_dRFP);
        GFPhcomp.AddConstant('k_dGFP', k_dGFP);

        GFPhcomp.AddConstant('K_R1',K_R1);
        GFPhcomp.AddConstant('K_R2',K_R2);
        GFPhcomp.AddConstant('K_R1R3',K_R1R3);
        GFPhcomp.AddConstant('K_R3R1',K_R3R1);
        GFPhcomp.AddConstant('K_R2R4',K_R2R4);
        GFPhcomp.AddConstant('K_R4R2',K_R4R2);
        GFPhcomp.AddConstant('K_R5',K_R5);

        GFPhcomp.AddConstant('n_R1', n_R1);
        GFPhcomp.AddConstant('n_R2', n_R2);
        GFPhcomp.AddConstant('n_R3', n_R3);
        GFPhcomp.AddConstant('n_R4', n_R4);
        GFPhcomp.AddConstant('n_R5', n_R5);

        dR1dt = GFPhcomp.AddCompositor('R1', input_species_1);
        dR2dt = GFPhcomp.AddCompositor('R2', input_species_2);
        dR3dt = GFPhcomp.AddCompositor('R3', 0);
        dR4dt = GFPhcomp.AddCompositor('R4', 0);
        dR5dt = GFPhcomp.AddCompositor('R5', 0);
        dBFPdt = GFPhcomp.AddCompositor('BFP', 0);
        dRFPdt = GFPhcomp.AddCompositor('RFP', 0);
        dGFPdt = GFPhcomp.AddCompositor('GFP', 0);

        GFPhcomp.AddPart(Part('NOT1', dR4dt,...
            Rate('k_R4*( ((K_R1)^(n_R1))/( ((K_R1)^(n_R1))+((R1)^(n_R1)) ) ) - k_dR4*R4')));

        GFPhcomp.AddPart(Part('NOT2', dR3dt,...
            Rate('k_R3*( ((K_R2)^(n_R2))/( ((K_R2)^(n_R2))+((R2)^(n_R2)) ) ) - k_dR3*R3')));

        GFPhcomp.AddPart(Part('NOR1', [dR5dt, dBFPdt],...
            [Rate('k_R5*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) )'), ...
            Rate('k_BFP*( (((K_R1R3)^(n_R1))/( ((K_R1R3)^(n_R1))+((R1)^(n_R1)) ))*(((K_R3R1)^(n_R3))/( ((K_R3R1)^(n_R3))+((R3)^(n_R3))) ) ) - k_dBFP*BFP')]));

        GFPhcomp.AddPart(Part('NOR2', [dR5dt, dRFPdt],...
            [Rate('k_R5*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dR5*R5'), ...
            Rate('k_RFP*( (((K_R2R4)^(n_R2))/( ((K_R2R4)^(n_R2))+((R2)^(n_R2)) ))*(((K_R4R2)^(n_R4))/( ((K_R4R2)^(n_R4))+((R4)^(n_R4))) ) ) - k_dRFP*RFP')]));

        GFPhcomp.AddPart(Part('NOT3', dGFPdt,...
            Rate('k_GFP*(((K_R5^n_R5)/((K_R5^n_R5)+(R5^n_R5)))^2) - k_dGFP*GFP')));
        
        [T,Y] = GFPhcomp.run([0,500]);
        GFP_output_high(n,m)=Y(end,GFPhcomp.CompositorIndex('GFP'));
    end
end

GFP_threshold_high = 296.8;
GFP_threshold_low = 299.2;
RFP_threshold = 298.9;
BFP_threshold = 298.9;

for i = 1:1:length(BFP_output(:,1))
    for k = 1:1:length(BFP_output(1,:))
        if BFP_output(i,k) >= BFP_threshold
            BFP_output(i,k) = 1;
        elseif BFP_output(i,k) < BFP_threshold
            BFP_output(i,k) = 0;
        end
        
        if RFP_output(i,k) >= RFP_threshold
            RFP_output(i,k) = 1;
        elseif RFP_output(i,k) < RFP_threshold
            RFP_output(i,k) = 0;
        end
 
        if GFP_output_high(i,k) >= GFP_threshold_high
            GFP_output_high(i,k) = 1;
        elseif GFP_output_high(i,k) < GFP_threshold_high
            GFP_output_high(i,k) = 0;
        end
        
        if GFP_output_low(i,k) >= GFP_threshold_low
            GFP_output_low(i,k) = 1;
        elseif GFP_output_low(i,k) < GFP_threshold_low
            GFP_output_low(i,k) = 0;
        end
    end
end

figure()
ax1 = subplot(2,2,1);
contourf(x_scaling_BFP,y_scaling_BFP,BFP_output)
map5 = [1, 1, 1
    0, 0, 1];
colormap(ax1, map5)
%colorbar
%xlabel('R1 molecules/cell'); 
ylabel('R2 molecules/cell'); 
%title('BFP Output (mlcl/cell)')

ax2 = subplot(2,2,4);
contourf(x_scaling_RFP,y_scaling_RFP,RFP_output)
map6 = [1, 1, 1
    1, 0, 0];
colormap(ax2, map6)
%colorbar
xlabel('R1 molecules/cell'); 
%ylabel('R2 molecules/cell'); title('RFP Output (mlcl/cell)')

ax3 = subplot(2,2,3);
contourf(x_scaling_GFPl,y_scaling_GFPl,GFP_output_low)
map7 = [1, 1, 1
    0, 1, 0];
colormap(ax3, map7)
%colorbar
xlabel('R1 molecules/cell'); 
ylabel('R2 molecules/cell'); 
%title('GFP Output (mlcl/cell)')

ax4 = subplot(2,2,2);
contourf(x_scaling_GFPh,y_scaling_GFPh,GFP_output_high)
colormap(ax4, map7)
%colorbar
%xlabel('R1 molecules/cell'); ylabel('R2 molecules/cell'); title('GFP Output (mlcl/cell)')
