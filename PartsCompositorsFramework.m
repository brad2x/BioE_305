%%PARTS AND COMPOSITORS FRAMEWORK
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