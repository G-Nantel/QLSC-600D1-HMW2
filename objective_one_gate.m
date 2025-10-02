function err = objective_one_gate(p, t, vStep, iData, Er, edgeMask)
% Parameters: [log10(gBar_mS)  Vhalf_m  b_m  tau_m]
p = enforce_one_gate(p);
[iModel, ~] = simulate_one_gate(p, t, vStep, Er);
res = (iModel - iData);
res(~edgeMask,:) = 0;                   % ignore step-edge transients
err = sum(res(:).^2);                   % SSE across all sweeps
end

function p = enforce_one_gate(p)
% Make sure physically required params are positive where needed
p(1) = max(-8, min(2, p(1)));           % log10(gBar) bounded
p(3) = max(1e-3, p(3));                 % slope >= 0 (increasing with V)
p(4) = max(0.1, p(4));                  % tau >= 0.1 ms
end
