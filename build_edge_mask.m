function keep = build_edge_mask(t, vStep, window_ms)
% Returns a logical mask (same size as vStep) that is true away from step edges
% so the optimizer ignores capacitive spikes at ~step on/off.
dt = t(2)-t(1);
w = max(1, round(window_ms/dt));   % half-window in samples
[nT, nS] = size(vStep);
keep = true(nT, nS);
for s = 1:nS
    v = vStep(:,s);
    edges = find(abs(diff(v)) > 2);           % detect step transitions
    idx = [];
    for e = edges(:)'
        idx = [idx, max(1, e-w):min(nT, e+w)]; %#ok<AGROW>
    end
    keep(unique(idx), s) = false;
end
end
