clc;
clear;
close all;


f = @(t,y) (-50*y) + sin(t); 


t0 = 0;
tf= 2;
y0 = 1;      
h = 0.04;
N = (tf - t0) / h;
      

t = linspace(t0, tf, N+1);


y_forward = zeros(1, N+1);
y_modified = zeros(1, N+1);
y_backward = zeros(1, N+1);
y_rk2 = zeros(1, N+1);
y_adams_bashforth = zeros(1, N+1);


y_forward(1) = y0;
y_modified(1) = y0;
y_backward(1) = y0;
y_rk2(1) = y0;
y_adams_bashforth(1) = y0;


% Forward Euler Method
for i = 1:N
    y_forward(i+1) = y_forward(i) + h * f(t(i), y_forward(i));
end

% Heun's Method
for i = 1:N
    y_pred = y_modified(i) + h * f(t(i), y_modified(i)); 
    y_modified(i+1) = y_modified(i) + (h / 2) * (f(t(i), y_modified(i)) + f(t(i+1), y_pred)); 
end

% Backward Euler Method 
for i = 1:N
    y_backward(i+1) = (y_backward(i) + h * sin(t(i+1))) / (1 + 50*h); 
end

% Runge-Kutta 2nd Order Method 
for i = 1:N
    k1 = h * f(t(i), y_rk2(i));
    k2 = h * f(t(i) + h, y_rk2(i) + k1);
    y_rk2(i+1) = y_rk2(i) + (k1 + k2) / 2;
end

% Adams-Bashforth Method 
y_adams_bashforth(2) = y_adams_bashforth(1) + h * f(t(1), y_adams_bashforth(1));
for i = 2:N
    y_adams_bashforth(i+1) = y_adams_bashforth(i) + (h/2) * (3*f(t(i), y_adams_bashforth(i)) - f(t(i-1), y_adams_bashforth(i-1)));
end


% Plot Results
figure;
plot(t, y_forward, 'r', 'LineWidth', 1.5); hold on;
plot(t, y_modified, 'b', 'LineWidth', 1.5);
plot(t, y_backward, 'g', 'LineWidth', 1.5);
plot(t, y_rk2, 'm', 'LineWidth', 1.5); 
plot(t, y_adams_bashforth, 'c', 'LineWidth', 1.5);
legend('Forward Euler', 'Modified Euler', 'Backward Euler','Runge-Kutta 2nd Order', 'Adams-Bashforth');
xlabel('Time t'); ylabel('y(t)');
title('Methods');
grid on;
