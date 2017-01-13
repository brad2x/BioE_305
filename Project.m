urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
         'part_composition_setup.m');
rehash;
part_composition_setup('v5');

learning = BioSystem();

%Dilution and Degradation Constants
learning.AddConstant('k_dil',0);
learning.AddConstant('k_deg',2);

%Protein production terms (assuming fast RNA dynamics)
learning.AddConstant('kp_R1',100);
learning.AddConstant('kp_R2',100);
learning.AddConstant('kp_R3',100);
learning.AddConstant('kp_TetR',100);
learning.AddConstant('kp_R4',100);
learning.AddConstant('kp_RFP',100);
learning.AddConstant('kp_R5',100);
learning.AddConstant('kp_R6',100);
learning.AddConstant('kp_Out',100);

%Michaelis-Menten constants
learning.AddConstant('K_R1',10);
learning.AddConstant('K_R2',10);
learning.AddConstant('K_R4',10);
learning.AddConstant('K_TetR',10);
learning.AddConstant('K_R2',10);
learning.AddConstant('K_R5',10);
learning.AddConstant('K_R6',10);

%Sequestration rate constants
learning.AddConstant('k_on_S1',50);
learning.AddConstant('k_off_S1',5);
learning.AddConstant('k_on_S2',50);
learning.AddConstant('k_off_S2',5);
learning.AddConstant('k_on_TetR',50);
learning.AddConstant('k_off_TetR',5);

%Compositors::..
%Inputs: S1, S2, Reset (Dox)
%Outputs: RFP, Out
%Protein-only Intermediates: R3, R4, R5, R6
%Protein-SM Intermediates: R1, S1R1, R2, S2R2, TetR, TetRDox
dS1dt = learning.AddCompositor('S1',0);
dS2dt = learning.AddCompositor('S2',0);
dR1dt = learning.AddCompositor('R1',100);
dR2dt = learning.AddCompositor('R2',100);
dS1R1dt = learning.AddCompositor('S1R1',0);
dS2R2dt = learning.AddCompositor('S2R2',0);
dDoxdt = learning.AddCompositor('Dox',0);
dR3dt = learning.AddCompositor('R3',0);
dTetRdt = learning.AddCompositor('TetR',0);
dTetRDoxdt = learning.AddCompositor('TetRDox',0);
dR4dt = learning.AddCompositor('R4',1000);
dRFPdt = learning.AddCompositor('RFP',0);
dR5dt = learning.AddCompositor('R5',0);
dR6dt = learning.AddCompositor('R6',1000);
dOutdt = learning.AddCompositor('Out',0);

%Parts::
%Parts logic occurs at the gate level, only including necessary
%intermediates
%1) Input Sensing
%2) Conditioning Gate
%3) Toggle Memory
%4) Recall Gate
%5) Output Gate

%learning.AddPart(...
%    Part('Input Sensing',[dS1dt dR1dt dS1R1dt dR3dt dS2dt dR2dt dS2R2dt], ...
%    [ Rate('k_off_S1*S1R1 - k_on_S1*S1*R1'), Rate('kp_R1 + k_off_S1*S1R1 - k_on_S1*S1*R1 - k_deg*R1'),...
%    Rate('-k_off_S1*S1R1 + k_on_S1*S1*R1'),...
%    Rate('kp_R3*(K_R1/(K_R1 + R1))-k_deg*(R3)'),...
%    Rate('k_off_S2*S2R2 - k_on_S2*S2*R2'), Rate('kp_R2 + k_off_S2*S2R2 - k_on_S2*S2*R2 - k_deg*R2'),...
%    Rate('-k_off_S2*S2R2 + k_on_S2*S2*R2')] ) );

learning.AddPart(...
    Part('Input Sensing', [dS1dt dR1dt dR3dt dS2dt dR2dt], ...
    [ Rate('-k_on_S1'), Rate('kp_R1 - k_on_S1*S1 - k_deg*R1'),...
    Rate('kp_R3*(S1/(K_R1 + S1)) - k_deg*R3'), Rate('-k_on_S1'),...
    Rate('kp_R2 - k_on_S2*S2 - k_deg*R2') ] ) );

%[T,Y] = learning.run([0,1000]);

[T, Y] = learning.run_pulses([...
    Pulse(0, 'S1', 0), ... % initial conditions
    Pulse(250, 'S1', 100), ... % inducton
    Pulse(1000, 'S1', 0), ... % stop the simulation 
]);

figure();
plot(T, Y(:, learning.CompositorIndex('S1')), T, Y(:, learning.CompositorIndex('R1')),...
    T, Y(:, learning.CompositorIndex('S1R1')), T, Y(:, learning.CompositorIndex('R3')));...
    %T, Y(:, learning.CompositorIndex('S2')), T, Y(:, learning.CompositorIndex('R2')),...
    %T, Y(:, learning.CompositorIndex('S2R2')));


xlabel('Time'); ylabel('Concentration'); legend('S1','R1','S1R1','R3','S2','R2','S2R2');

%%
learning.AddPart(...
    Part('Conditioning Gate',dTetRdt, ...
    Rate('kp_TetR*((K_R1/(K_R1 + R1)*(K_R2/(K_R2 + R2)) - k_deg*TetR') ) );

learning.AddPart(...
    Part('Memory Toggle',[dTetRdt dDoxdt dTetRDoxdt dR4dt dRFPdt], ...
    [ Rate('kp_TetR*(K_R4/(K_R4 + R4)) + k_off_TetR*(TetRDox) - k_on_TetR*(TetR)*(Dox) - k_deg*TetR '),... 
    Rate('k_off_TetR*(TetRDox) - k_on_TetR*(TetR)*(Dox) - k_deg*Dox'),...
    Rate('k_on_TetR*(TetR)*(Dox) - k_off_TetR*(TetRDox)'),...
    Rate('kp_R4*(K_TetR/(K_TetR + TetR)) - k_deg*(R4)'),...
    Rate('kp_RFP*(K_R4/(K_R4 + R4) - k_deg*(RFP)') ] ) );

learning.AddPart(...
    Part('Recall Gate',dR5dt, ...
    Rate('kp_R5*((K_R4/(K_R4 + R4))*(K_R2/(K_R2 + R2)) - k_deg*(R5)') ) );

learning.AddPart(...
    Part('Output Gate',[dR6dt dOutdt], ...
    [ Rate('kp_R6*( (K_R5/(K_R5 + R5))*(K_R3/(K_R3 + R3)) ) - k_deg*(R6)'),...
    Rate('kp_Out*(K_R6/(K_R6 + R6)) - k_deg*(Out)') ] ) );


[T,Y] = learning.run([0,1000]);


