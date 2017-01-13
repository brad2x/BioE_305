%Problem Set 2

% this block fetches the framework: do not edit
urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
    'part_composition_setup.m');
rehash;
part_composition_setup('v5');
% the system is constructed & simulated in the block below: feel free to edit

% define a BioSystem: this object will hold the parts and compositors
sys = BioSystem();

% define constants. assume that concentration units are in nM and time in s

sys.AddConstant(Const('k_if', 1000));
sys.AddConstant(Const('k_ir', 100));
sys.AddConstant(Const('k_arf', 10));
sys.AddConstant(Const('k_arr', 1));
sys.AddConstant(Const('k_rf', 10));
sys.AddConstant(Const('k_rr', 1));
sys.AddConstant(Const('k_oc1', 0.1));
sys.AddConstant(Const('k_oc2', 0.1));
sys.AddConstant(Const('k_tln', 0.1));
sys.AddConstant(Const('k_mdeg', 0.1)); % mRNA degradation & dilution
sys.AddConstant(Const('k_pdeg', 0.01)); % protein degradation & dilution
sys.AddConstant(Const('k_f', 10)); % rate of activator binding DNA
sys.AddConstant(Const('k_r', 1)); % rate of activator unbinding DNA

%Compositors for Protein 1
dInddt = sys.AddCompositor('Ind', 0); % inducer, set value later
dactdt = sys.AddCompositor('act', 5); % nM, about 5 molecules in E. coli
dActdt = sys.AddCompositor('Act', 0);
dDNA1dt = sys.AddCompositor('DNA1',1);
dActDNA1dt = sys.AddCompositor('ActDNA1',0);
dRNAPdt = sys.AddCompositor('RNAP', 1); % nM, about 1 copy in an E. coli
dActDNA1RNAPdt = sys.AddCompositor('ActDNA1RNAP',0);
dmRNA1dt = sys.AddCompositor('mRNA1',0);
dP1dt = sys.AddCompositor('P1', 0);

%Additional Compositors for Protein 2
dDNA2dt = sys.AddCompositor('DNA2', 0.5); %0.5
dDNA2RNAPdt = sys.AddCompositor('DNA2RNAP',0);
dmRNA2dt = sys.AddCompositor('mRNA2',0);
dP2dt = sys.AddCompositor('P2', 0);
% TODO: define parts and add them to the system
% Remember to use the same parts as in the Groundwork section of this problem set

Protein1 = Part('Protein 1',...
    [dactdt dInddt dActdt dDNA1dt dActDNA1dt dRNAPdt dActDNA1RNAPdt dmRNA1dt dP1dt],...
    [Rate('-k_if*act*Ind + k_ir*Act') 
    Rate('-k_if*act*Ind + k_ir*Act') 
    Rate('k_if*act*Ind - k_ir*Act - k_f*Act*DNA1 + k_r*ActDNA1') 
    Rate('-k_f*Act*DNA1 + k_r*ActDNA1') 
    Rate('k_f*Act*DNA1 - k_r*ActDNA1 - k_arf*ActDNA1*RNAP + k_arr*ActDNA1RNAP + k_oc1*ActDNA1RNAP') 
    Rate('-k_arf*ActDNA1*RNAP + k_arr*ActDNA1RNAP + k_oc1*ActDNA1RNAP') 
    Rate('k_arf*ActDNA1*RNAP - k_arr*ActDNA1RNAP - k_oc1*ActDNA1RNAP') 
    Rate('k_oc1*ActDNA1RNAP - k_mdeg*mRNA1') 
    Rate('k_tln*mRNA1 - k_pdeg*P1')]);

Protein2 = Part('Protein 2', ...
[dDNA2dt dRNAPdt dDNA2RNAPdt dmRNA2dt dP2dt],...
[Rate('-k_rf*DNA2*RNAP + k_rr*DNA2RNAP + k_oc2*DNA2RNAP') 
    Rate('-k_rf*DNA2*RNAP + k_rr*DNA2RNAP + k_oc2*DNA2RNAP') 
    Rate('k_rf*DNA2*RNAP - k_rr*DNA2RNAP - k_oc2*DNA2RNAP') 
    Rate('k_oc2*DNA2RNAP - k_mdeg*mRNA2') 
    Rate('k_tln*mRNA2 - k_pdeg*P2')]);

sys.AddPart(Protein1);
sys.AddPart(Protein2);

% solve/simulate the system. change the amount of inducer at 3 times
[T, Y] = sys.run_pulses([...
    Pulse(0, 'Ind', 0), ...     % initial conditions
    Pulse(1000, 'Ind', 10), ... % spike in 10 nM of inducer @ t = 1000
    Pulse(2000, '', 0), ...     % stop the simulation at time 2000
]);


%plot
figure();
plot(T, Y(:, sys.CompositorIndex('P1')), ...
     T, Y(:, sys.CompositorIndex('P2')), ...
     T, Y(:, sys.CompositorIndex('RNAP')), ...
     T, Y(:, sys.CompositorIndex('Ind')))
ylim([0 12]);
legend('P1', 'P2', 'RNAP', 'Ind', 'Location', 'best')
xlabel('Time (s)');
ylabel('Concentration (nM)');

%figure();
%plot(T,Y(:, sys.CompositorIndex('P2')));
%xlabel('Time'); ylabel('P2');

ratio = (Y(132,13)-Y(end,13))/ Y(132,13);

% plot all compositors, useful for debugging:
% figure();
% num_compositors = size(Y, 2);
% for i = 1:num_compositors
%     subplot(num_compositors, 1, i);
%     plot(T, Y(:, i));
%     xlabel('Time');
%     ylabel(sprintf('[%s]', sys.compositors(i).name), 'Rotation', 90);
% end

%%
%%Positive Autoregulation
urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
    'part_composition_setup.m');
rehash;
part_composition_setup('v5');

sys = BioSystem();

sys.AddConstant(Const('k_if', 1000));

dInddt = sys.AddCompositor('Ind', 0); % inducer, set value later
