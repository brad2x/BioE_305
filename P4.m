%% TODO
%   1.  Fill in the parameter values from the problem specification
%       and the parts using the ODEs from the previous problem. This
%       problem is graded on the rates and parts being correct.
%   2.  Select a parameter or species initial value to iterate over
%       for the plotting. The problems below will ask about this.

        
%% This block fetches the framework: do not edit
urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
    'part_composition_setup.m');
rehash;
part_composition_setup('v5');


%% The system is constructed & simulated in the block below: feel free to edit
%   With the provided initial conditions and the parameters and rates described
%   in the problem description, you should get oscillations.

% Do not change the name
oscillator = BioSystem();


% You can edit this to increase/decrease the simulation run time
simulation_run_time = 2000;


% Species + initial conditions
A = oscillator.AddCompositor('A', 10); 
B = oscillator.AddCompositor('B', 0); 
m_A = oscillator.AddCompositor('m_A', 0);
m_B = oscillator.AddCompositor('m_B', 0);
p = oscillator.AddCompositor('p', 5);
pA = oscillator.AddCompositor('pA', 0);


% Constants 
%   TODO: Add these from problem description
%Recall you add constants in the form: oscillator.AddConstant('', NaN);
oscillator.AddConstant('K_A', 1);
oscillator.AddConstant('K_B', 2);
oscillator.AddConstant('n_A', 2);
oscillator.AddConstant('n_B', 4);
oscillator.AddConstant('k_tln', 5);
oscillator.AddConstant('k_txn', 5);
oscillator.AddConstant('k_mdeg', 0.5);
oscillator.AddConstant('k_pdeg', 0.05);
oscillator.AddConstant('k_on', 50);
oscillator.AddConstant('k_off', 50);

% Parts
%   TODO: Add these from problem description
oscillator.AddPart(Part('All interactions', [A B m_A m_B p pA], ...
    [Rate('k_tln * m_A - k_pdeg * A - k_on * n_A * A^n_A * p + k_off * n_A * pA'), ...
     Rate('k_tln*m_B-k_pdeg*B'), ...
     Rate('k_txn*((A^n_A)/((A^n_A)+(K_A^n_A)))*((K_B^n_B)/(( K_B^n_B)+(B^n_B)))-k_mdeg*m_A'), ...
     Rate('k_txn*((A^n_A)/((A^n_A)+(K_A^n_A)))-k_mdeg*m_B'), ...
     Rate('-k_on*A^n_A*p+k_off*pA'), ...
     Rate('k_on*A^n_A*p-k_off*pA')]));

%% Plotting iteration
%   TODO: Indicate which variable or parameter you would
%   like to change with each loop, plotting the changes.
%   
%   To start, we recommend looking at m_B and iterating p.
%
%   You should only iterate over one species or one 
%   parameter at a time. If you do not want to iterate
%   over one or the other, set it to false.

plot_species = 'm_B';       % The species to be plotted
iter_species = false;         % The species to be iterated with different initial values
iter_constant = 'k_off';      % The parameter to be iterated with different values
iter_values = 50:25:200;          % The values to iterate over


%% Plotting/simulation, we recommend not changing the below
%   Simulate system from 0 to simulation_run_time units

figure()
ax = gca();
hold(ax, 'on')
ax.ColorOrder = jet(length(iter_values));
for val = iter_values
    
    if (any(iter_species))
        oscillator.ChangeInitialValue(iter_species, val);
        iterable = iter_species;
    elseif (any(iter_constant))
        oscillator.ChangeConstantValue(iter_constant, val);
        iterable = iter_constant;
    else
        % Do not change any initial conditions or paramters
        iterable = '';
    end
    
    [t, y] = oscillator.run([0 simulation_run_time]);

    % Plot concentration of species [speciesIndex] over time
    plot(ax, t, y(:, oscillator.CompositorIndex(plot_species)))
end
xlabel('Time (s)');
ylabel(sprintf('Concentration of %s', plot_species))
split_iterable = strsplit(iterable, '_');
if (length(split_iterable) == 2)
    % Contains an underscore
    iterable = strcat(split_iterable{1}, '_{', split_iterable{2}, '}');
end
legend(arrayfun(@(x) sprintf('%s = %d', iterable, x), iter_values, 'unif', false), 'Location', 'Best')
title('An oscillator with load', 'fontsize', 14)
      