function r1gaussianMixtureFittingRaRb
    S.x = -5:.1:5;
    S.fh = figure('units','normalized',...
    'Position', [0.515 0.025 0.415 0.87],... %%%%
              'name','slider_plot');
     % Declar the member variables of class S
     S.ma = 0.5;
     S.mb = 1;
     S.ra = 3;
     S.rb = 3;
     
     S.maSlider = uicontrol('style','slider',...
                                       'unit','normalized',...
                                       'position',[0.2 0.1 0.7 0.01],...
                                       'min',0,'max',5,...
                                       'sliderstep',[0.1, 0.1],...
                                       'callback', {@gaussianMixtureValue, 'ma'}); 
                                   
     S.mbSlider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 0.15 0.7 0.01],...
               'min',0,'max',3,...
               'sliderstep',[0.1, 0.1],...
               'callback', {@gaussianMixtureValue, 'mb'});
           
     S.raSlider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 0.2 0.7 0.01],...
               'min',0,'max',10,...
               'sliderstep',[1, 1],...
               'callback', {@gaussianMixtureValue, 'ra'});
           
     S.rbSlider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 0.25 0.7 0.01],...
               'min',0,'max',10,...
               'sliderstep',[1, 1],...
               'callback', {@gaussianMixtureValue, 'rb'});
    
     guidata(S.fh, S);  

end

function gaussianMixtureValue(hSlider, eventData, param)
S = guidata(hSlider); 
if(param == 'ma')
%     disp("ma value is ")
%     get(hSlider, 'value')
    S.ma = get(hSlider, 'value');
elseif(param == 'mb')
%     disp("mb value is ")
%     get(hSlider, 'value')
    S.mb = get(hSlider, 'value'); 
elseif(param == 'ra')
    S.ra = get(hSlider, 'value'); 
elseif(param == 'rb')
    S.rb = get(hSlider, 'value'); 
end
x = -5:0.01:5;

z = rectangularPulse(-1, 1, x);
w = 1/exp(-0.5) *  exp(-x.^2/2);
y = S.ra * S.ma * exp(-S.mb * x.^2/2)/sqrt(2*pi) - S.ma * exp(-S.rb * S.mb * x.^2/2)/sqrt(2*pi);

plot(x, y, 'g')
hold on
plot(x, z, 'b')
plot(x, w, 'c')
hold off

disp("ra value is ")
S.ra
disp("rb value is ")
S.rb
disp("mb value is ")
S.mb
disp("ma value is ")
S.ma

guidata(S.fh, S);
end
