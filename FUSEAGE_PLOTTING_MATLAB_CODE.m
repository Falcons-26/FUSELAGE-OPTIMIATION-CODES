function find_k_critical_from_merged_csv()
clc;

TARGET_STRESS = 7.58/3.0;

fprintf("=== Find k where MY/I = 7.58/3 (using single merged CSV + airfoil) ===\n");

%% --- Input merged CSV ---
merged_file = input("Enter merged MOI+M CSV file path [default results_moi_with_M.csv]: ","s");
if isempty(merged_file)
    merged_file = "results_moi_with_M.csv";
end
if ~isfile(merged_file)
    error("File not found: %s", merged_file);
end

%% --- Input airfoil ---
airfoil_file = input("Enter airfoil profile CSV/TXT file path [default airfoil.csv]: ","s");
if isempty(airfoil_file)
    airfoil_file = "airfoil.csv";
end
if ~isfile(airfoil_file)
    error("File not found: %s", airfoil_file);
end

%% --- Read merged CSV ---
df = readtable(merged_file);

required = ["x","k","centroid","I","M"];

if ~ismember("h_eff", df.Properties.VariableNames)
    if all(ismember(["y_total","baseplate_line"], df.Properties.VariableNames))
        df.h_eff = df.y_total - df.baseplate_line;
    else
        error("Merged CSV must include h_eff OR (y_total + baseplate_line)");
    end
end
required = [required,"h_eff"];

for i = 1:length(required)
    df.(required(i)) = double(df.(required(i)));
end

df = df(~any(ismissing(df(:,required)),2),:);
df = df(abs(df.I) > 1e-12,:);
if isempty(df)
    error("No valid rows after cleaning.");
end

%% --- Read airfoil ---
[~,~,ext] = fileparts(airfoil_file);
ext = lower(ext);

if ext == ".txt"
    airfoil = readmatrix(airfoil_file);
    airfoil = array2table(airfoil(:,1:2),"VariableNames",["x","y"]);
else
    airfoil = readtable(airfoil_file);
end

airfoil.x = double(airfoil.x);
airfoil.y = double(airfoil.y);
airfoil = airfoil(~any(ismissing(airfoil),2),:);

%% --- Split upper / lower surfaces ---
x_vals = airfoil.x;
last_dec = height(airfoil);
for i = 2:length(x_vals)
    if x_vals(i) >= x_vals(i-1)
        last_dec = i;
        break
    end
end

upper    = airfoil(1:last_dec-1,:);
lower_af = airfoil(last_dec:end,:);   % ✅ renamed

if upper.x(1) > upper.x(end), upper = flipud(upper); end
if lower_af.x(1) > lower_af.x(end), lower_af = flipud(lower_af); end

[xu,yu] = collapse_duplicates(upper.x, upper.y);
[xl,yl] = collapse_duplicates(lower_af.x, lower_af.y);

y_top_at    = @(xq) interp_safe(xq,xu,yu);
y_bottom_at = @(xq) interp_safe(xq,xl,yl);

%% --- Find k-critical ---
results = table();

unique_x = unique(df.x);
for xi = 1:length(unique_x)
    x_val = unique_x(xi);
    group = sortrows(df(df.x==x_val,:), "k");

    ytop = y_top_at(x_val);
    ybot = y_bottom_at(x_val);
    thickness_at_x = ytop - ybot;

    mask = group.k >= 0 & group.k <= group.h_eff/3.5;
    if ~isnan(thickness_at_x)
        mask = mask & group.k <= thickness_at_x;
    end
    group = group(mask,:);
    if isempty(group), continue; end

    c_top = group.h_eff - group.centroid;
    group = group(c_top > 0,:);
    c_top = c_top(c_top > 0);

    stress = (group.M .* c_top) ./ group.I;
    [~,idx] = min(abs(stress - TARGET_STRESS));
    best = group(idx,:);

    y_neutral = ytop - best.k;

    results = [results; table( ...
        x_val, best.k, y_neutral, stress(idx), best.M, ...
        best.centroid, best.I, abs(stress(idx)-TARGET_STRESS), ...
        'VariableNames',{'x','k_critical','y_neutral','stress','M','centroid','I','abs_error'})];
end

results = sortrows(results,"x");

%% --- Exclusion range ---
exclude_range = input("Enter exclusion x range as 'start,end' (blank for none): ","s");
exclude_start = NaN; exclude_end = NaN;
if ~isempty(exclude_range)
    v = str2double(strsplit(exclude_range,","));
    if numel(v)==2
        exclude_start = min(v);
        exclude_end   = max(v);
    end
end

df_plot = results;
if ~isnan(exclude_start)
    df_plot = df_plot(~(df_plot.x>=exclude_start & df_plot.x<=exclude_end),:);
end

%% --- Baseplate & offset ---
if ismember("baseplate_line", df.Properties.VariableNames)
    baseplate_max = max(df.baseplate_line,[],'omitnan');
else
    xs = linspace(min(airfoil.x),max(airfoil.x),200);
    baseplate_max = max(arrayfun(y_bottom_at,xs),[],'omitnan');
end

x_offset_val = input("Enter x offset above baseplate (mm) [default 10]: ");
if isempty(x_offset_val), x_offset_val = 10; end
x_horizontal_y = baseplate_max + x_offset_val;

%% --- Save ---
out_file = input("Enter output CSV filename [default results_k_critical.csv]: ","s");
if isempty(out_file)
    out_file = "results_k_critical.csv";
end
writetable(results,out_file);
fprintf("✅ Saved %s\n", out_file);

%% --- Plot ---
figure; hold on; grid on;
plot(airfoil.x, airfoil.y,'Color',[0.6 0.6 0.6],'LineWidth',1.5);
plot(df_plot.x, df_plot.y_neutral,'ro-','LineWidth',1.2);
plot(results.x, arrayfun(y_top_at,results.x),'g--');
plot(results.x, arrayfun(y_bottom_at,results.x),'b--');
yline(baseplate_max,'k','LineWidth',1.5);
yline(x_horizontal_y,'m--','LineWidth',1.2);
xlabel("x (mm)");
ylabel("y (mm)");
title("Airfoil outline with neutral axis vs x");
legend("Airfoil","Neutral axis","Top","Bottom","Baseplate","x-offset");
axis tight;

end

%% ---------- helpers ----------
function [xp2,fp2] = collapse_duplicates(xp,fp)
keep = [true; diff(xp)>0];
xp2 = xp(keep);
fp2 = fp(keep);
end

function y = interp_safe(xq,xp,fp)
if numel(xp)<2
    y = NaN;
else
    y = interp1(xp,fp,xq,'linear','extrap');
end
end
