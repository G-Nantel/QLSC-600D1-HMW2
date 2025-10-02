% Cook Assignment 1 — minimal ONE-GATE model (constant tau)
% ----------------------------------------------------------
% Loads data, fits 4 params, and makes the required 5‑panel figure.

clear; clc; close all;

% Data: t [ms], vStep [mV], iUnknownCurrent [mA]
load('CookAssignemnt1UnknownCurrent.mat');
Er = -15;                                                               % mV (given)  :contentReference[oaicite:3]{index=3}
dt = t(2)-t(1);

% --- Tiny helper: ignore ±2 ms around step edges (capacitive spikes)
nT = numel(t); nS = size(vStep,2);
keep = true(nT,nS);
edges = [zeros(1,nS); abs(diff(vStep))] > 2;                            % detect jumps
w = max(1, round(2/dt));                                                % 2 ms half-window
for s=1:nS, ei = find(edges(:,s));
    for e = ei', keep(max(1,e-w):min(nT,e+w),s) = false; end
end

% --- Parameters (one gate): p = [log10(gBar_mS), Vhalf_mV, log(b), log(tau_ms)]
p0 = [log10(5e-4), -35, log(0.08), log(10)];                            % sensible seeds  :contentReference[oaicite:4]{index=4}
obj = @(p) sum(sum(((simulate_one_gate(p,t,vStep,Er) - iUnknownCurrent).*keep).^2));
p  = fminsearch(obj, p0, optimset('Display','off'));                    % Nelder–Mead

% --- Final simulation with fitted params
[iModel, x, g, xInf] = simulate_one_gate(p, t, vStep, Er);

% --- Required figure ----------------------------------------------------
fh = figure('Color','w','Position',[100 100 900 1000]);
sgtitle('Cook H&H-like current: one-gate fit');

% (1) x_inf(V)
subplot(5,1,1);
V = linspace(-80,20,400);
plot(V, xInf(V), 'LineWidth', 2); grid on; ylim([-0.05 1.05]);
xlabel('mV'); ylabel('steady-state');
txt = sprintf('gBar = %.6g mS\nEr = %g mV\nV1/2 = %.3g mV, b = %.3g, tau = %.3g ms', ...
              10.^p(1), Er, p(2), exp(p(3)), exp(p(4)));
annotation('textbox',[0.62 0.78 0.3 0.16],'String',txt,'FitBoxToText','on','EdgeColor','none');

% (2) voltage steps
subplot(5,1,2); plot(t, vStep, 'LineWidth', 1.1); grid on; ylabel('V step (mV)');

% (3) gate over time
subplot(5,1,3); plot(t, x, 'LineWidth', 1.1); grid on; ylabel('x_a (activation)');

% (4) conductance
subplot(5,1,4); plot(t, g, 'LineWidth', 1.1); grid on; ylabel('g (mS)');

% (5) currents: data (solid) vs model (dashed)
subplot(5,1,5);
plot(t, iUnknownCurrent, 'LineWidth', 2); hold on;
set(gca,'ColorOrderIndex',1);
plot(t, iModel, '--', 'LineWidth', 2); grid on; xlabel('msec'); ylabel('i (mA)');
legend('Data','Model','Location','SouthEast');

% --- Export single‑page PDF (as requested in the handout)
teamName = 'INFINITE_LOOP';                                                  % <--- change me
print(fh, sprintf('%s_CookHHLab1.pdf', teamName), '-dpdf','-painters'); % :contentReference[oaicite:5]{index=5}

% =========================== Local function =============================
function [i, x, g, xInf] = simulate_one_gate(p, t, vStep, Er)
    % p = [log10(gBar), Vhalf, log(b), log(tau)]  (b>0, tau>0 by construction)
    gBar = 10.^p(1); Vhalf = p(2); b = exp(p(3)); tau = exp(p(4));
    dt = t(2)-t(1); alpha = 1 - exp(-dt/tau);

    xInf = @(V) 1 ./ (1 + exp(-b*(V - Vhalf)));                         % logistic  :contentReference[oaicite:6]{index=6}
    nT = numel(t); nS = size(vStep,2);
    x = zeros(nT,nS); g = zeros(nT,nS); i = zeros(nT,nS);

    for s=1:nS
        v = vStep(:,s);
        x(1,s) = xInf(v(1));                                            % start at steady-state
        for k=2:nT
            x(k,s) = x(k-1,s) + (xInf(v(k)) - x(k-1,s))*alpha;          % X(t+dt)=X+...(1-exp(-dt/tau))
        end
        g(:,s) = gBar * x(:,s);
        i(:,s) = g(:,s) .* (v - Er);
    end
end
