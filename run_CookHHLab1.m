clear; clc; close all;

% Load Data
load('CookAssignemnt1UnknownCurrent.mat');   % t [ms], vStep [mV], iUnknownCurrent [mA]
Er = -15;                                     % reversal potential (mV) (given)
dt = t(2) - t(1);                             % sample interval (ms)

% Edge-spike masking (capacitive transients)
nT = numel(t); nS = size(vStep,2);
keep = true(nT,nS);
edges = [false(1,nS); abs(diff(vStep)) > 2];  % detect step transitions
w = max(1, round(2/dt));                      % ±2 ms half-window
for s = 1:nS
    ei = find(edges(:,s));
    for e = ei'
        keep(max(1,e-w):min(nT,e+w), s) = false;  % ignore edge samples
    end
end

% ----------------------- Parameterization --------------------------------
% p = [log10(gBar_mS), Vhalf_mV, b_raw, tau_raw]
% Map b, tau smoothly into realistic bounds so the fit stays physical.
b_min = 0.03;  b_max = 0.15;      % mV^-1
t_min = 5;     t_max = 30;        % ms
sigm  = @(z) 1./(1+exp(-z));

p0  = [log10(5e-4), -35, 0, 0];   % seeds: b,tau at mid-bounds

obj = @(p) sum(((simulate_one_gate(p,t,vStep,Er,b_min,b_max,t_min,t_max) ...
                 - iUnknownCurrent).*keep).^2, 'all');

p   = fminsearch(obj, p0, optimset('Display','off'));

[iModel, x, g, xInf, par_out] = simulate_one_gate(p, t, vStep, Er, ...
                                                  b_min, b_max, t_min, t_max);
b_disp   = par_out.b;           % mapped values for display
tau_disp = par_out.tau;
gBar_disp= par_out.gBar;
Vhalf    = par_out.Vhalf;

% Figure
fh = figure('Color','w','Position',[100 80 980 1040]);
sgtitle('Cook H&H-like current: one-gate fit (masked edges, bounded params)');

% (1) x_inf(V)
subplot(5,1,1);
V = linspace(-80,20,400);
plot(V, xInf(V), 'LineWidth', 2); grid on; ylim([-0.05 1.05]);
xlabel('mV'); ylabel('steady-state');

% Parameter textbox only (no equations)
txt = sprintf('gBar = %.6g mS\nEr = %g mV\nV_{1/2} = %.3g mV,  b = %.3g mV^{-1},  \\tau = %.3g ms', ...
              gBar_disp, Er, Vhalf, b_disp, tau_disp);
annotation('textbox',[0.62 0.78 0.32 0.16],'String',txt, ...
           'FitBoxToText','on','EdgeColor','none','Interpreter','tex');

% (2) voltage steps
subplot(5,1,2);
plot(t, vStep, 'LineWidth', 1.1); grid on; ylabel('V step (mV)');

% (3) gate over time
subplot(5,1,3);
plot(t, x, 'LineWidth', 1.1); grid on; ylabel('x_a (activation)');

% (4) conductance
subplot(5,1,4);
plot(t, g, 'LineWidth', 1.1); grid on; ylabel('g (mS)');

% (5) currents: Data (solid) vs Model (dashed)
subplot(5,1,5);
plot(t, iUnknownCurrent, 'LineWidth', 2); hold on;
set(gca,'ColorOrderIndex',1);
plot(t, iModel, '--', 'LineWidth', 2);
grid on; xlabel('msec'); ylabel('i (mA)');
legend('Data','Model','Location','SouthEast');

% ----------------------------- Export ------------------------------------
teamName = 'INFINITE_LOOP';   % <-- change to your team name
print(fh, sprintf('%s_CookHHLab1.pdf', teamName), '-dpdf','-painters');

% =========================== Local function ==============================
function [i, x, g, xInf, pars] = simulate_one_gate(p, t, vStep, Er, ...
                                                   b_min, b_max, t_min, t_max)
    % Map parameters -> physical ranges
    sigm = @(z) 1./(1+exp(-z));
    gBar  = 10.^p(1);                          % mS
    Vhalf = p(2);                               % mV
    b     = b_min + (b_max - b_min) * sigm(p(3));  % mV^-1
    tau   = t_min + (t_max - t_min) * sigm(p(4));  % ms

    dt    = t(2) - t(1);
    alpha = 1 - exp(-dt / tau);                % discrete-time easing factor

    xInf = @(V) 1 ./ (1 + exp(-b * (V - Vhalf)));  % steady-state sigmoid

    [nT, nS] = deal(numel(t), size(vStep,2));
    x = zeros(nT,nS); g = zeros(nT,nS); i = zeros(nT,nS);

    for s = 1:nS
        v = vStep(:,s);
        x(1,s) = xInf(v(1));                   % start at steady-state
        for k = 2:nT
            x(k,s) = x(k-1,s) + (xInf(v(k)) - x(k-1,s)) * alpha;
            % x(k,s) = min(1, max(0, x(k,s))); % optional safety clamp
        end
        g(:,s) = gBar * x(:,s);
        i(:,s) = g(:,s) .* (v - Er);
    end

    pars = struct('gBar', gBar, 'Vhalf', Vhalf, 'b', b, 'tau', tau);
end
