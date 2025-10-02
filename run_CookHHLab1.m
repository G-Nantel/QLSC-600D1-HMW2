% QLSC 600 — Cook Assignment 1
% Driver: loads data, fits model, produces required figure & PDF
% -------------------------------------------------------------------------
clear; clc; close all;

% --- Load data from the assignment (t, vStep, iUnknownCurrent)
%     File variables are described in the PDF handout.
load('CookAssignemnt1UnknownCurrent.mat');  % t [ms], vStep [mV], iUnknownCurrent [mA]

Er = -15;           % mV  (given in the assignment "Hints") 
modelVariant = 1;   % 1 = one-gate m;

% --- Build a mask to de-emphasize capacitive spikes around step edges
edgeMask = build_edge_mask(t, vStep, 2.0);  % exclude +/- 2 ms around detected step edges

% --- Initial parameter guesses (sensible based on the slides/example)
% One-gate: [log10(gBar_mS)  Vhalf_m   slope_b   tau_m_ms]
p0_1 = [log10(5e-4),  -35,  0.08,  10];

% Two-gate: [log10(gBar_mS) Vhalf_m bm tau_m  Vhalf_h bh tau_h]
p0_2 = [log10(6e-4),  -35,  0.08,  6,   -55, -0.06,  20];

% --- Fit
opts = optimset('Display','iter','MaxFunEvals',2e4,'MaxIter',2e4);
switch modelVariant
    case 1
        obj = @(p) objective_one_gate(p, t, vStep, iUnknownCurrent, Er, edgeMask);
        pFit = fminsearch(obj, p0_1, opts);
    case 2
        obj = @(p) objective_two_gate(p, t, vStep, iUnknownCurrent, Er, edgeMask);
        pFit = fminsearch(obj, p0_2, opts);
end

% --- Simulate with the fitted parameters
switch modelVariant
    case 1
        [iModel, comp, pars] = simulate_one_gate(pFit, t, vStep, Er);
    case 2
        [iModel, comp, pars] = simulate_two_gate(pFit, t, vStep, Er);
end

% --- Make the required figure (matches the example layout in the handout)
fh = figure('Color','w','Position',[100 100 900 1000]);

% (a) m∞(V) / h∞(V)
subplot(5,1,1);
Vgrid = linspace(-80, 20, 400);
plot(Vgrid, comp.mInfHandle(Vgrid), 'LineWidth', 2); hold on;
if isfield(comp,'hInfHandle')
    plot(Vgrid, comp.hInfHandle(Vgrid), 'LineWidth', 2);
    legend('m_\infty','h_\infty','Location','NorthWest');
else
    legend('x_\infty (m)','Location','NorthWest');
end
xlabel('mV'); ylabel('steady-state'); ylim([-0.05 1.05]); grid on;

% Parameter textbox
txt = compose_params_text(pars, Er);
annotation('textbox',[0.62 0.78 0.3 0.16],'String',txt,'FitBoxToText','on','EdgeColor','none');

% (b) V step
subplot(5,1,2);
plot(t, vStep, 'LineWidth', 1.5);
ylabel('V step (mV)'); grid on; xlim([t(1) t(end)]);

% (c) Gates over time (per sweep)
subplot(5,1,3);
plot(t, comp.m, 'LineWidth', 1.5); hold on;
if isfield(comp,'h')
    plot(t, comp.h, '--', 'LineWidth', 1.5);
    ylabel('m  (solid),  h  (dashed)');
else
    ylabel('x_a (activation)');
end
grid on; xlim([t(1) t(end)]);

% (d) g(t)
subplot(5,1,4);
plot(t, comp.g, 'LineWidth', 1.5);
ylabel('g (mS)'); grid on; xlim([t(1) t(end)]);

% (e) Current: data (solid) vs model (dashed)
subplot(5,1,5);
plot(t, iUnknownCurrent, 'LineWidth', 2); hold on;
set(gca,'ColorOrderIndex',1);           % repaint same colors for model
plot(t, iModel, '--', 'LineWidth', 2);
xlabel('msec'); ylabel('i (mA)'); grid on; xlim([t(1) t(end)]);
legend('Data','Model','Location','SouthEast');

sgtitle('Cook H&H-like current: data vs. fitted model');

% --- Export a single-page PDF you can submit
teamName = 'INFINITE-LOOP';   % <-- change this
outpdf = sprintf('%s_CookHHLab1.pdf', teamName);
set(fh, 'PaperPositionMode', 'auto');
print(fh, outpdf, '-dpdf', '-painters');

fprintf('\nSaved figure to: %s\n', outpdf);
