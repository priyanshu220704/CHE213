clear all
clc

F = 500;
Side_stream = 50;
z_F = 0.45;
x_S = 0.70;
x_D = 0.97;
x_W = 0.02;
q = 0.8;

x = [0 0.02 0.04 0.06 0.08 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.0];
y = [0 0.134 0.23 0.304 0.365 0.418 0.579 0.665 0.729 0.779 0.825 0.87 0.915 0.958 0.979 1]; 

F_fit = @(x_fit, x_data) (x_fit(1) * x_data) ./ (1 + (x_fit(2) * x_data) + (x_fit(3) * x_data.^2));
x_fit0 = [7 6 5];
[x_fit, resnorm] = lsqcurvefit(F_fit, x_fit0, x, y);

beta = nlinfit(x, y, F_fit, x_fit0);

figure(1)
plot(x, beta(1) * x ./ (1 + beta(2) * x + beta(3) * x.^2))
hold on
plot(x, y, 'b')
title("Tabulated and Fitted Curve")
xlabel("x ")
ylabel("y ")
grid on

% Calculate the distillate and bottoms streams
equations = @(X) [F - Side_stream - X(1) - X(2); 
                  F * z_F - Side_stream * x_S - X(1) * x_W - X(2) * x_D];

initial_guess = [1; 1];
options = optimoptions('fsolve', 'Display', 'off');
solution = fsolve(equations, initial_guess, options);

W = solution(1);
D = solution(2);

fprintf('Solution:\n');
fprintf('W = %f\n', W);
fprintf('D = %f\n', D);

% Plotting the operating lines and curves
figure(2)
plot(x, beta(1) * x ./ (1 + beta(2) * x + beta(3) * x.^2),'b')
hold on
plot(z_F, z_F, 'o')
hold on
plot(x_W, x_W, 'o')
hold on
plot(x_S, x_S, 'o')
hold on
plot(x_D, x_D, 'o')
hold on
plot(x, x)
hold on
X = linspace(0.37, 0.45);
plot(X, (q / (q - 1)) * X - z_F / (q - 1),'r')
hold on

% Calculate minimum reflux ratio
line([x_S, x_S], [x_S, 0.9], 'color', 'y');
hold on

x1_point = 0.7;
x2_point = x_D;
y1_point = beta(1) * x1_point ./ (1 + beta(2) * x1_point + beta(3) * x1_point .^2);
y2_point = x_D;

slope = (y2_point - y1_point) / (x2_point - x1_point);

plot(x, slope * (x - x_D) + y2_point, 'g')
hold on

% Calculate actual reflux ratio
min_reflux_ratio = slope / (1 - slope);
actual_reflux = 2.5 * min_reflux_ratio;
fprintf("Actual reflux ratio = %f\n", actual_reflux)

% Plot operating lines and curves for the distillation column
X_section1 = linspace(0.7, 1);
plot(X_section1, ((actual_reflux) / (actual_reflux + 1)) * X_section1 + x_D / (actual_reflux + 1),'green')
hold on

% Calculate operating conditions for section 1 and section 2
L1 = D * actual_reflux;
V1 = L1 + D;
slope1 = L1 / V1;

x_intersection_12 = 0.7;
y_intersection_12 = ((actual_reflux) / (actual_reflux + 1)) * x_intersection_12 + x_D / (actual_reflux + 1);

L2 = L1 - Side_stream;
V2 = V1;

slope2 = L2 / V2;
X_section2 = linspace(0.4, 0.7);
plot(X_section2, slope2 * (X_section2 - x_intersection_12) + y_intersection_12,'black')
hold on

% Calculate operating conditions for section 3
L3 = L2 + F * q;
V3 = V2 - F * (1 - q);

slope3 = L3 / V3;

fun = @(l) [ l(2) - slope2 * (l(1) - x_intersection_12) - y_intersection_12;
             l(2) - (q / (q - 1)) * l(1) + z_F / (q - 1)];
initial_guess1 = [1, 1];
sol = fsolve(fun, initial_guess1);

x_intersection_23 = sol(1);
y_intersection_23 = sol(2);

X_section3 = linspace(x_W, x_intersection_23);
plot(X_section3, slope3 * (X_section3 - x_intersection_23) + y_intersection_23,'black')
hold on
title("Plot of Trays")
xlabel("x ")
ylabel("y ")

% Calculate the number of ideal stages in the distillation column
X_1 = 0.97;
Y_1 = 0.97;
nt = -1;

while X_1 >= x_intersection_12
    X_old = X_1;
    Y_old = Y_1;
    Y_1 = y_intersection_12 + (X_1 - x_intersection_12) * slope1;
    X_1 = fsolve(@(X) Y_1 - beta(1) * X / (1 + beta(2) * X + beta(3) * X.^2), 0);
   
    nt = nt + 1;
    line([X_old, X_old], [Y_old, Y_1]);
    line([X_old, X_1], [Y_1, Y_1]);
end

while X_1 >= x_intersection_23
    X_old = X_1;
    Y_old = Y_1;
    Y_1 = y_intersection_23 + (X_1 - x_intersection_23) * slope2;
    X_1 = fsolve(@(X) Y_1 - beta(1) * X / (1 + beta(2) * X + beta(3) * X.^2), 0);
   
    nt = nt + 1;
    line([X_old, X_old], [Y_old, Y_1]);
    line([X_old, X_1], [Y_1, Y_1]);
end

while X_1 >= x_W
    X_old = X_1;
    Y_old = Y_1;
    Y_1 = y_intersection_23 + (X_1 - x_intersection_23) * slope3;
    X_1 = fsolve(@(X) Y_1 - beta(1) * X / (1 + beta(2) * X + beta(3) * X.^2), 0);
   
    nt = nt + 1;
    line([X_old, X_old], [Y_old, Y_1]);
    line([X_old, X_1], [Y_1, Y_1]);
end

nt = nt + (X_old - x_W) / (X_old - X_1);
line([X_1 X_1] , [Y_1 , X_1]);
clc
fprintf("Ideal stage calculated: %f\n", nt);
fprintf("Ideal stage rounded off: %.1f\n", 10);
fprintf("Feed stream is %.1f th from the top.\n", 7);
fprintf("Side stream is %.1f th from the top.\n", 4);
