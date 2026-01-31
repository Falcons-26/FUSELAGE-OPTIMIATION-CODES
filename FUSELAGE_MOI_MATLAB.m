function fuselage_moi_merged_runner()
clc;

%% ------------------ Macaulay bending moment ------------------
fprintf("=== Fuselage MOI + Merged M Generator (MATLAB) ===\n");

%% --- Read airfoil CSV/TXT ---
csv_file = input("Enter CSV/TXT file path (first column x, second column y): ","s");
if ~isfile(csv_file)
    error("File not found.");
end

try
    raw = readmatrix(csv_file);
    raw = raw(:,1:2);
    raw = raw(~any(isnan(raw),2),:);
    x_raw = raw(:,1);
    y_raw = raw(:,2);
catch
    error("Could not parse airfoil file.");
end

if isempty(x_raw)
    error("Airfoil file empty after cleaning.");
end

%% --- Baseplate gap ---
baseplate_gap = input("Enter minimum vertical gap above lower surface (mm) [default 8]: ");
if isempty(baseplate_gap)
    baseplate_gap = 8;
end

%% --- Geometry inputs ---
w     = input("Enter flange width w (mm): ");
tw    = input("Enter flange thickness tw (mm): ");
rt    = input("Enter web thickness rt (mm): ");
x_bot = input("Enter bottom web height x_bot (mm): ");
k_min = input("Enter minimum k (mm): ");
step  = input("Enter step size for k sweep (mm): ");

%% --- Point loads ---
n_loads = input("Enter number of point loads: ");
point_loads = zeros(n_loads,2);

for i = 1:n_loads
    point_loads(i,1) = input(sprintf("Enter magnitude of load %d (N): ",i));
    point_loads(i,2) = input(sprintf("Enter position of load %d (mm): ",i));
end

%% --- Split airfoil into upper / lower ---
forward_end = length(x_raw);
for i = 2:length(x_raw)
    if x_raw(i) >= x_raw(i-1)
        forward_end = i;
        break
    end
end

upper = [x_raw(1:forward_end-1), y_raw(1:forward_end-1)];
lower = [x_raw(forward_end:end), y_raw(forward_end:end)];

% Ensure increasing x
if upper(1,1) > upper(end,1), upper = flipud(upper); end
if lower(1,1) > lower(end,1), lower = flipud(lower); end

[xu,yu] = collapse_duplicates(upper(:,1), upper(:,2));
[xl,yl] = collapse_duplicates(lower(:,1), lower(:,2));

fixed_x = min(x_raw);

results = [];

fprintf("\nDEBUG: x | y_top | y_bot | thickness | baseplate | h_eff\n");

%% --- Main loop ---
for i = 1:length(xu)
    x_val = xu(i);
    y_top = yu(i);

    % nearest neighbor bottom
    [~,idx] = min(abs(xl - x_val));
    y_bot = yl(idx);

    thickness = abs(y_top) + abs(y_bot);
    baseplate_line = y_bot + baseplate_gap;
    h_eff = y_top - baseplate_line;

    fprintf("%.3f | %.3f | %.3f | %.3f | %.3f | %.3f\n", ...
            x_val, y_top, y_bot, thickness, baseplate_line, h_eff);

    if h_eff <= 0
        continue
    end

    % bending moment
    dist = max(0, x_val - fixed_x);
    M = macaulay_bending_moment(point_loads, dist);

    local_k_max = h_eff - (tw + x_bot);
    if local_k_max < k_min
        continue
    end

    [k_list, centroids, inertias] = ...
        sweep_centroid_inertia(h_eff, w, tw, rt, x_bot, k_min, step, local_k_max);

    for j = 1:length(k_list)
        results = [results; ...
            x_val, y_top, y_bot, thickness, baseplate_line, h_eff, ...
            k_list(j), centroids(j), inertias(j), M];
    end
end

%% --- Save CSV ---
out_file = input("Enter output CSV filename [default results_moi_with_M.csv]: ","s");
if isempty(out_file)
    out_file = "results_moi_with_M.csv";
end

headers = {'x','y_top','y_bottom','y_total','baseplate_line','h_eff','k','centroid','I','M'};
T = array2table(results,'VariableNames',headers);
writetable(T,out_file);

fprintf("\nAll results saved to %s\n", out_file);
fprintf("Baseplate gap = %.3f mm\n", baseplate_gap);
fprintf("Done.\n");

end

%% ================= HELPER FUNCTIONS =================

function M = macaulay_bending_moment(loads, x)
M = 0;
for i = 1:size(loads,1)
    if x >= loads(i,2)
        M = M + loads(i,1)*(x - loads(i,2));
    end
end
end

function [k_list, centroid_list, inertia_list] = ...
    sweep_centroid_inertia(h, w, tw, rt, x_bot, k_min, step, k_max_override)

k_upper = min(k_max_override, h/2);
if k_upper < k_min
    k_list = []; centroid_list = []; inertia_list = [];
    return
end

A_f = w * tw;
y_f = tw / 2;

k_list = [];
centroid_list = [];
inertia_list = [];

k = k_min;
while k <= k_upper + 1e-9

    A_b = 2 * rt * x_bot;
    A_t = 2 * rt * k;
    A_tot = A_f + A_b + A_t;

    y_b = tw + x_bot/2;
    y_t = h - k/2;

    ybar = (A_f*y_f + A_b*y_b + A_t*y_t) / A_tot;

    I_f = (w*tw^3)/12 + A_f*(ybar - y_f)^2;
    I_s = (rt/6)*(x_bot^3 + k^3) ...
          + 2*rt*(x_bot*(y_b-ybar)^2 + k*(y_t-ybar)^2);

    I_tot = I_f + I_s;

    k_list = [k_list; k];
    centroid_list = [centroid_list; ybar];
    inertia_list = [inertia_list; I_tot];

    k = k + step;
end
end

function [xp2,fp2] = collapse_duplicates(xp,fp)
keep = [true; diff(xp) > 0];
xp2 = xp(keep);
fp2 = fp(keep);
end
