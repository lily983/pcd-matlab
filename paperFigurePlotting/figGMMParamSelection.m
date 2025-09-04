close all; clc;clear

[a1, b1]= get_gmm_param(eye(1), '1SG');

[a2, b2]=get_gmm_param(eye(1), '5SG');

x = -3:0.01:4;

psi_2gm =a1(1) * exp(-b1(1) * x.^2/2)/sqrt(2*pi) + a1(2)* exp(-b1(2) * x.^2/2)/sqrt(2*pi);

psi_5gm =a2(1) * exp(-b2(1) * x.^2/2)/sqrt(2*pi) + a2(2)* exp(-b2(2) * x.^2/2)/sqrt(2*pi) + a2(3)* exp(-b2(3) * x.^2/2)/sqrt(2*pi)...
    +a2(4)* exp(-b2(4) * x.^2/2)/sqrt(2*pi) + a2(5)* exp(-b2(5) * x.^2/2)/sqrt(2*pi);

inta = rectangularPulse(-1, 1, x);

% Single Gaussian function f(x)
k = sqrt(2*pi)/exp(-0.5);
fx = k*exp(-x.^2/2)/sqrt(2*pi);


%% set y limit plot
figure; hold on;
Linewidth = 2.5;
% plot(x, z, 'color', hex2rgb('D95980'), 'LineWidth', Linewidth);
plot(x, psi_2gm, 'color', hex2rgb('284E60'), 'LineWidth', Linewidth);
plot(x, psi_5gm, 'color', hex2rgb('F99B45'), 'LineWidth', Linewidth);
plot(x, fx, 'color', hex2rgb('63AAC0'), 'LineWidth', Linewidth);
% set(gca, 'YScale', 'log')
ylim([0.0, 10])
xlim([-3 3]);

%% zoom plot
figure; hold on;
Linewidth = 2;
plot(x, inta, 'color', hex2rgb('D95980'), 'LineWidth', Linewidth);
plot(x, psi_2gm, 'color', hex2rgb('284E60'), 'LineWidth', Linewidth);
plot(x, psi_5gm, 'color', hex2rgb('F99B45'), 'LineWidth', Linewidth);
plot(x, fx, 'color', hex2rgb('63AAC0'), 'LineWidth', Linewidth);
ylim([1, 1e+09])
xlim([-2 4]);
set(gca, 'YScale', 'log')

%%
% Get current axis position and limits
p = gca;

% set zoomin figure position
pos=[0.1, -0.3, 0.4, 0.4];

% Calculate x,y points of zoomPlot
x1 = (pos(1)-p.Position(1))/p.Position(3)*diff(p.XLim)+p.XLim(1);
x2 = (pos(1)+pos(3)-p.Position(1))/p.Position(3)*diff(p.XLim)+(p.XLim(1));
y1 = (pos(2)-p.Position(2))/p.Position(4)*diff(p.YLim)+p.YLim(1);
y2 = ((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(p.YLim)+p.YLim(1);

% set xbounds=[x1 x2] zoom indicies
xbounds = [-2, 2];

% Plot lines connecting zoomPlot to original plot points
index = find(x>=xbounds(1) & x<=xbounds(2)); % Find indexes of points in zoomPlot
% rectangle('Position',[xbounds(1) min(psi_2gm(index)) diff(xbounds) max(psi_2gm(index))-min(psi_2gm(index))]);
hold on
plot([xbounds(1) x1], [0 y1], 'k')
plot([xbounds(2) x2], [0 y1], 'k')

% Plot zoomPlot and change axis
z = axes('position',pos);
box on 
hold on
z.YLim = [0,2];
z.XLim = xbounds;
% plot(x,y)
Linewidth2 = 1.5;
plot(x, inta, 'color', hex2rgb('D95980'), 'LineWidth', Linewidth2);
plot(x, psi_5gm, 'color', hex2rgb('F99B45'), 'LineWidth', Linewidth2);
plot(x, fx, 'color', hex2rgb('63AAC0'), 'LineWidth', Linewidth2);
plot(x, psi_2gm, 'color', hex2rgb('284E60'), 'LineWidth', Linewidth2);

% adjust plot parameters
p.FontName = 'Times New Roman';
z.FontName = 'Times New Roman';

p.FontSize = 15
z.FontSize = 15

% p.XAxisLocation = 'origin';
% p.YAxisLocation = 'origin';
%%
figure;hold on;box on 
% plot(x,y)
Linewidth2 = 1.5;
plot(x, inta, 'color', hex2rgb('D95980'), 'LineWidth', Linewidth2);
plot(x, psi_5gm, 'color', hex2rgb('F99B45'), 'LineWidth', Linewidth2);
plot(x, fx, 'color', hex2rgb('63AAC0'), 'LineWidth', Linewidth2);
plot(x, psi_2gm, 'color', hex2rgb('284E60'), 'LineWidth', Linewidth2);
ylim([0,2])
xlim([-2,2])