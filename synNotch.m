urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
    'part_composition_setup.m');
rehash;
part_composition_setup('v5');

synSys = BioSystem();

sysSyn.AddConstant(Const('Kd_sn-a',100));
sysSyn.AddConstant(Const('Vsn-a', 100);
sysSyn.AddConstant(Const('Vcleave', 100);
sysSyn.AddConstant(Const('Vmrna1', 100));
sysSyn.AddConstant(Const('Vd_tf-re', 100);
sysSyn.AddConstant(Const('Kmrna1-deg', 100);


dSynNotchdt = sysSyn.AddCompositor('SynNotch',0);
dAntigendt = sysSyn.AddCompositor('Antigen',0);
dSynNotch:Antigendt = sysSyn.AddCompositor('SynNotch:Antigen',0);
dSynNotchTFdt = sysSyn.AddCompositor('SynNotchTF',0);
dSynNotchREdt = sysSyn.AddCompositor('SynNotchRE', 0);
dmRNA1dt= sysSyn.AddCompositor('mRNA1', 0);


sysSyn.AddPart('SynNotch + Antigen -K_sna-> SynNotch:Antigen',...
    [dSynNotchdt dAntigendt dSynNotch:Antigendt],...
    [rate('(-1*Vsn-a*Antigen)/(Kd_sn-a + Antigen)'), rate('-1*Vsn-a*Antigen)/(Kd_sn-a + Antigen)'),...
    rate('(Vsn-a*Antigen)/(Kd_sn-a + Antigen)')]);
sysSyn.AddPart('SynNotch:Antigen -Kcleave-> SynNotchTF', ...
    [dSynNotch:Antigendt dSynNotchTFdt], ...
    [rate('-1*Vcleave*SynNotch:Antigen') rate('Vcleave*SynNotch:Antigen')]);
sysSyn.AddPart('SynNotchTF -ktf-deg-> 0', [dSynNotchTFdt], []);
sysSyn.AddPart('SynNotchTF + SynNotchRE --> mRNA1', [dSynNotchTFdt dSynNotchREdt dmRNA1dt],...
    [rate('0'), rate('0'), rate('(Vmrna1*SynNotchTF)/(Kd_tf-re + SynNotchTF)')]);
sysSynAddPart('mRNA1 -Kmrna1-deg-> 0', [dmRNA1dt], [rate('-kmrna1-deg*mRNA1')]);


[T,Y] = sys.run([0,25]);

figure();
plot(T, Y(:, sys.CompositorIndex('X')));
xlabel('Time'); ylabel('Concentration of X');

