% A Bit of Practice: Plotting with MATLAB
figure();           %Create a new figure
cells_in_culture = @(t,N) (log(2)/30)*N; %Expression for dN/dt as an 
% "anonymous" function where t = time, N = concentration.

tspan = 0:5:120;

hold on;

[T_1, N_1] = ode23s(cells_in_culture, tspan, 1);
[T_2, N_2] = ode23s(cells_in_culture, tspan, 2);
[T_10, N_10] = ode23s(cells_in_culture, tspan, 10);
plot(T_1, N_1, 'b',T_2,N_2,'g',T_10, N_10, 'r');
xlabel('Time (minutes)'); ylabel('Cells per unit volume');
legend('1 cell','2 cells','10 cells');

%% Number of Cells in a Finite Culture: Equilibria
K = 300; 
k = log(2)/30;

dNdt = @(t,N) k*N*(1-N/K);

[t_0, N_0] = ode23s(dNdt, [0,600], 0);
[t_1, N_1] = ode23s(dNdt, [0,600], 1);
[t_10, N_10] = ode23s(dNdt, [0,600], 10);
[t_300, N_300] = ode23s(dNdt, [0,600], 300);
[t_400, N_400] = ode23s(dNdt, [0,600], 400);

figure(); hold on;
plot(t_0,N_0,t_1, N_1, t_10, N_10, t_300, N_300, t_400, N_400);
xlabel('Time (minutes)'); ylabel('Cells per unit volume');
legend('nothing','1 cell','10 cells','300 cells','400 cells');


%% Parts and Compositors Framework
urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
    'part_composition_setup.m');
rehash;
part_composition_setup('v5');

sys = BioSystem();
% Create and add constants and compositors
sys.AddConstant('k', 0.1);
dAdt = sys.AddCompositor('A',10);   %State variable representing
% rate of change of [A], which starts at A(0) = 10
dBdt = sys.AddCompositor('B',0);
dEdt = sys.AddCompositor('E',1);

%Define and Add the parts
reaction = Part('A + E -k> B + E', [dAdt dBdt dEdt], [Rate('-k*A*E') Rate('k*A*E') Rate('0')]);
% this process involves A, B, and E
% how this process affects A, B, E
sys.AddPart(reaction);

%simulate the system from t = 0 to t = 25
[T,Y] = sys.run([0,25]);
%T holds a vector of time and Y is a matrix where each row is a vector of
%the values of the system variables corresponding to the compositors (in
%order of addition) for each moment of time in T

%plot the amount of B vs. time
figure();
plot(T, Y(:, sys.CompositorIndex('B')));
xlabel('Time'); ylabel('Concentration of B');