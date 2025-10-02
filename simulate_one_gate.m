function [iModel, comp, pars] = simulate_one_gate(p, t, vStep, Er)
% p = [log10(gBar_mS), Vhalf_m, b_m, tau_m]
gBar = 10.^p(1);
Vhalf = p(2);
bm = p(3);
tau_m = p(4);

dt = t(2) - t(1);
alpha = 1 - exp(-dt / tau_m);

mInf = @(V) 1 ./ (1 + exp(-bm*(V - Vhalf)));  % logistic

nT = numel(t);
nS = size(vStep,2);
m = zeros(nT, nS);
g = zeros(nT, nS);
iModel = zeros(nT, nS);

for s = 1:nS
    v = vStep(:,s);
    m(1,s) = mInf(v(1));                % start at steady-state for initial V
    for k = 2:nT
        m(k,s) = m(k-1,s) + (mInf(v(k)) - m(k-1,s))*alpha;
    end
    g(:,s) = gBar * m(:,s);
    iModel(:,s) = g(:,s) .* (v - Er);
end

% components for plotting/report
comp.m = m;
comp.g = g;
comp.mInfHandle = mInf;
pars = struct('gBar_mS', gBar, 'Vhalf_m', Vhalf, 'b_m', bm, 'tau_m_ms', tau_m);
end
