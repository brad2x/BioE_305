%% TODO
%   1.  Write the rate equations for the decay of each species
%       Then run the simulation to demonstrate you have correctly
%       accounted for degradation/dilution for each species
%   2.  Look at the problems below and perturb the system to
%       answer them. You can run the MATLAB code either in this
%       window (with the RUN CODE button), in your own MATLAB 
%       environment, or in the MATLAB sandbox. 


%% Do not edit this:
urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
         'part_composition_setup.m');
rehash;
part_composition_setup('v5');


%% Feel free to edit the values of these constants to answer the questions below
% Constants were chosen based on realistic biological values
% Do not change the names of the constants or you will not get
% the correct answer for this problem itself

% Questions focus on these constants
k_dil = log(2) / 22;            % dilution due to cell division (per hour)
k_deg = log(2) / 46;             % protein degradation rate (per hour)
input_induction_time = 24;       % Time after T0 when I0 is added to system (hr)

% We suggest not changing these
k_trx = 3;                      % (pM per hour)
k_tln = 300;                    % (pM per hour)
k_mdeg = 9;                     % (per hour)
k_p = k_trx * k_tln / k_mdeg;   % maximum production rate (pM per hour)
k_basal = 0.001 * k_p;          % basal expression rate (pM per hour)
n = 2.0;                        % Hill coefficient for all binding interactions (AU)
K_0 = 1e3;                      % binding constant for inducer molecule (I0) (pM)
K_Act = 10e3;                   % binding constant for activator (Act) (pM)


%% Do not edit the following:

% The simulation is run with enough time for y = A * (1 - e^(-kt)) to reach 99.9% of A (max)
simulation_run_time = log(0.001) * -1 / min(k_deg, k_dil);
induction_value = 100 * K_0;

% Biosystem setup
transfection = BioSystem();
transfection.AddConstant('k_dil', k_dil);
transfection.AddConstant('k_basal', k_basal);
transfection.AddConstant('k_p', k_p);
transfection.AddConstant('n', n);
transfection.AddConstant('k_deg', k_deg);
transfection.AddConstant('K_0', K_0);
transfection.AddConstant('K_Act', K_Act);

dIOdt = transfection.AddCompositor('I0', 0);
dActdt = transfection.AddCompositor('Act', 0);
dMarkerdt = transfection.AddCompositor('Marker', 0);
dOutdt = transfection.AddCompositor('Out', 0);
dPlasmiddt = transfection.AddCompositor('Plasmid', 100);

transfection.AddPart(...
    Part('I0->Act', [dIOdt dActdt], [Rate('0'), ...
        Rate('k_basal * Plasmid + 2 * k_p * Plasmid * I0^n / (K_0^n + I0^n)')]));
transfection.AddPart(...
    Part('Act->Out', [dActdt dOutdt], [Rate('0'), ...
        Rate('k_basal * Plasmid + k_p * Plasmid * Act^n / (K_Act^n + Act^n)')]));
transfection.AddPart(...
    Part('->Marker', dMarkerdt, Rate('k_basal * Plasmid + 5 * k_p * Plasmid')));


%% Add parts for decay of Act, Marker, Out, and Plasmid
%   You must account for dilution and degradation accordingly
%   Only use the constants k_deg and k_dil, do not make any new ones
%       or your answer will be graded as wrong.
%   Assume 1st order dilution and degradation rates
%   Every species dilutes with rate k_dil

% Act degrades at rate k_deg
transfection.AddPart(...
    Part('Act->', dActdt, Rate('-1*k_deg*Act - k_dil*Act')));

% Marker degrades at rate 2 * k_deg
transfection.AddPart(...
    Part('Marker->', dMarkerdt, Rate('-2*k_deg*Marker - k_dil*Marker')));

% Out degrades at rate k_deg
transfection.AddPart(...
    Part('Out->', dOutdt, Rate('-1*k_deg*Out - k_dil*Out')));

% Plasmids do not degrade
transfection.AddPart(...
    Part('Plasmid->', dPlasmiddt, Rate('-1*k_dil*Plasmid')));


%% Do not edit anything below

[T_dil, Y_dil] = transfection.run_pulses([...
    Pulse(0, 'I0', 0), ... % initial conditions
    Pulse(input_induction_time, 'I0', induction_value), ... % inducton
    Pulse(simulation_run_time, '', 0), ... % stop the simulation 
]);

transfection.ChangeConstantValue('k_dil', 0);

[T_no_dil, Y_no_dil] = transfection.run_pulses([...
    Pulse(0, 'I0', 0), ... % initial conditions
    Pulse(input_induction_time, 'I0', induction_value), ... % inducton
    Pulse(simulation_run_time, '', 0), ... % stop the simulation 
]);

close all

figure();
subplot(2, 1, 1);
plot(T_dil, Y_dil(:, transfection.CompositorIndex('Act')))
hold on
plot(T_dil, Y_dil(:, transfection.CompositorIndex('Marker')))
plot(T_dil, Y_dil(:, transfection.CompositorIndex('Out')))
plot(T_dil, Y_dil(:, transfection.CompositorIndex('Plasmid')))
plot(T_dil, Y_dil(:, transfection.CompositorIndex('I0')))
title('Dilution due to cell division')
xlabel('Time')
ylabel('Concentration (pM)')
legend('Act', 'Marker', 'Out', 'Plasmid', 'I0', 'Location', 'best')

subplot(2, 1, 2);
plot(T_no_dil, Y_no_dil(:, transfection.CompositorIndex('Act')))
hold on
plot(T_no_dil, Y_no_dil(:, transfection.CompositorIndex('Marker')))
plot(T_no_dil, Y_no_dil(:, transfection.CompositorIndex('Out')))
plot(T_no_dil, Y_no_dil(:, transfection.CompositorIndex('Plasmid')))
plot(T_no_dil, Y_no_dil(:, transfection.CompositorIndex('I0')))
title('No dilution due to cell division')
xlabel('Time')
ylabel('Concentration (pM)')
legend('Act', 'Marker', 'Out', 'Plasmid', 'I0', 'Location', 'best')

figure();
subplot(2, 1, 1);
marker = Y_dil(:, transfection.CompositorIndex('Marker'));
plot(T_dil, Y_dil(:, transfection.CompositorIndex('Act')) ./ marker, '-')
hold on
plot(T_dil, marker ./ marker, '-')
out_ratio_dil = Y_dil(:, transfection.CompositorIndex('Out')) ./ marker;
plot(T_dil, out_ratio_dil, '-')
title('Dilution due to cell division: normalized to Marker')
xlabel('Time')
ylabel('Concentration/Marker')
legend('Act/Marker', 'Marker/Marker', 'Out/Marker', 'Location', 'best')

subplot(2, 1, 2);
marker = Y_no_dil(:, transfection.CompositorIndex('Marker'));
plot(T_no_dil, Y_no_dil(:, transfection.CompositorIndex('Act')) ./ marker, '-')
hold on
plot(T_no_dil, marker ./ marker, '-')
out_ratio_no_dil = Y_no_dil(:, transfection.CompositorIndex('Out')) ./ marker;
plot(T_no_dil, out_ratio_no_dil, '-')
title('No dilution due to cell division: normalized to Marker')
xlabel('Time')
ylabel('Concentration/Marker')
legend('Act/Marker', 'Marker/Marker', 'Out/Marker', 'Location', 'best')

figure();
subplot(2, 1, 1);
plot(T_dil, out_ratio_dil / out_ratio_dil(transfection.time_to_index(T_dil, input_induction_time)))
hold on
plot(T_no_dil, out_ratio_no_dil / out_ratio_no_dil(transfection.time_to_index(T_no_dil, input_induction_time)))
title('On/off ratio of Out/Marker')
xlabel('Time')
ylabel('On/off ratio of Out/Marker')
legend('Dilution', 'No dilution', 'Location', 'best')
subplot(2, 1, 2);
% run_pulses() calls an ODE solver that automatically determines at what times
% to determine the solution at. This means that the two simulated systems will
% not have matching T vectors. We must interpolate the data in order to be able
% to calculate the difference between the two Out time traces
% element-by-element. In this case, the no dilution simulation has fewer
% datapoints than the simulation with dilution.
[T1, Y1, T2, Y2] = transfection.interpolate_traces(...
    T_dil, out_ratio_dil / out_ratio_dil(transfection.time_to_index(T_dil, input_induction_time)), ...
    T_no_dil, out_ratio_no_dil / out_ratio_no_dil(transfection.time_to_index(T_no_dil, input_induction_time)));
plot(T1, 1 - abs(Y2 - Y1) ./ Y2, '-')
hold on
plot([0 simulation_run_time], [0.95 0.95], '-')
title('Agreement of normalized on/off (relative to no dilution)')
xlabel('Time')
ylabel('Agreement')
legend('Agreement', '95% threshold', 'Location', 'best')
      