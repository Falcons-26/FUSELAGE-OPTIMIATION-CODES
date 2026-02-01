clc;
clear;

fprintf('=== Critical Section Thickness Calculator ===\n');

% Input permissible stress
stress = round(input('Enter permissible stress: '), 4);

fprintf('\n--- Section with Highest Applied Force ---\n');

% Inputs
F  = round(input('Enter maximum point load applied (F): '), 4);
a2 = round(input('Enter cut height (2a): '), 4);
b2 = round(input('Enter cut width (2b = q): '), 4);
r  = round(input('Enter corner radius (r, can be 0): '), 4);
p  = round(input('Enter p (maximum distance at cut midpoint): '), 4);

% ---------------- KT Calculation ----------------
if a2 == 0
    error('Cut width (2a) cannot be zero');
end

ratio_r = round(r / a2, 4);   % r / (2a)

% a and b doubled as requested
a = a2;
b = b2;
ratio_ab = round(a / b, 4);

C1 = round(14.815 - 22.308*sqrt(ratio_r) + 16.298*ratio_r, 4);
C2 = round(-11.201 - 13.789*sqrt(ratio_r) + 19.200*ratio_r, 4);
C3 = round(0.2020 + 54.620*sqrt(ratio_r) - 54.748*ratio_r, 4);
C4 = round(3.232 - 32.530*sqrt(ratio_r) + 30.964*ratio_r, 4);

KT = round( ...
      C1 ...
    + round(C2 * ratio_ab, 4) ...
    + round(C3 * ratio_ab^2, 4) ...
    + round(C4 * ratio_ab^3, 4), 4);

fprintf('\nStress concentration factor Kt: %.4f\n', KT);

% ---------------- Thickness Calculation ----------------
if p <= b2
    error('p must be greater than q');
end

t = round((KT * F) / ((p - b2) * stress), 4);
fprintf('Required thickness: %.4f\n', t);

% Nearest integer (engineering rounding)
final_thickness = floor(t + 0.5);
fprintf('Final thickness (nearest integer): %d\n', final_thickness);