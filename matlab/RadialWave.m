%% Cleaning
clc;
clear variables;
close all;
warning('off');

%% Parameters
meanR = 0;
stdR2 = 0.005;
g = 9.81;
h0 = 10;
c = sqrt(g*h0);

%% Initial perturbation
g = @(r) 0.3*exp(-(1/2)*((r-meanR).^2)/(stdR2));

%% Loop to get the expression of u(r, t) -> D(R1, i)
R = 0:0.01:0.75;
t = 0:0.001:0.1;

for R1 = 1:length(R)
    RCurrent = R(R1);    

    for i = 1:length(t)
        T = t(i);
        
        % Poisson representation (see PDE lecture 8)
        F = @(r, theta) ...
            g(sqrt((RCurrent+r.*cos(theta)).^2 + (r.*sin(theta)).^2))...
            ./(sqrt((c^2*T^2)./(r.^2)-1));
        
        I(i) = (1/(2*pi*c))*integral2(F, 0, c*T, 0, 2*pi);
        
        if(i ~= 1)
            D(R1, i) = (I(i) - I(i-1))/(t(i)-t(i-1));
        end
    end
end

%% Animation of the height radially
figure('Name', 'Animation of the height radially');
for i=1:length(t)
    plot(R, h0+D(:, i), '-b');
    textTime = text(0.5, 10.1, ['$t=$ ', num2str(t(i)), ' [s]'], ...
        'Color', 'red', 'FontSize', 20);
    set(textTime, 'interpreter', 'latex');
    set(gca, 'Units', 'normalized', ...
        'FontUnits', 'points',...
        'FontWeight','normal',...
        'FontSize', 15,...
        'FontName','Times');
    xlabel('$r$', 'interpreter', 'latex');
    ylabel('$h(r)$', 'interpreter', 'latex');
    grid on;
    
    xlim([0, max(R)]);
    ylim(h0+[-0.15, 0.3]);
    pause(0.1);
end

%% Animation of the height spatially
figure('Name', 'Animation of the height spatially');
% Convert the radial component to a 2D map
[X,Y] = meshgrid(0:0.01:1,0:0.01:1);
Rradial = sqrt((X-0.5).^2 + (Y-0.5).^2);
[m, n] = size(X);
for i = 1:m
    for j = 1:n
        for k = 1:length(t)
            index = find(R >= Rradial(i, j), 1);
            if(isnan(D(index, k)))
                Z(i, j, k) = 0;
            else
                Z(i, j, k) = D(index, k);
            end
        end
    end
end

% Display the result
for i = 1:length(t)
    s = surf(X, Y, h0 + Z(:, :, i));
    % view(90, 0);
    s.EdgeColor = 'none';
    colorbar;
    caxis(h0+[-0.1, 0.3]);
    set(gca, 'Units', 'normalized', ...
        'FontUnits', 'points',...
        'FontWeight','normal',...
        'FontSize', 15,...
        'FontName','Times');
    xlabel('$x$', 'interpreter', 'latex');
    ylabel('$y$', 'interpreter', 'latex');
    zlabel('$h(x, y)$', 'interpreter', 'latex');
    grid on;
    
    xlim([0, 1]);
    ylim([0, 1]);
    zlim(h0+[-0.1, 0.3]);
    
    pause(0.01);
end