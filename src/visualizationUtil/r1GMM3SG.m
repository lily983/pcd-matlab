function r1GMM3SG
    S.x = -5:.1:5;
    S.fh = figure('units','normalized',...
    'Position', [0.515 0.025 0.415 1.0],... %%%%
              'name','slider_plot');
          
     distance = 0.03;
    
     %Initial value
     a=[-9.5616 7.9123  4.4867];
     b=[7.1483 3.4981 2.6616];
     
     S.a1 = a(1); S.a2 = a(2);  S.a3 = a(3);  
     S.b1 = b(1); S.b2 = b(2); S.b3 = b(3); 

     S.a1Slider = uicontrol('style','slider',...
                                       'unit','normalized',...
                                       'position',[0.2 7*distance 0.7 0.01],...
                                       'min',-10,'max',20,...
                                       'sliderstep',[0.1, 0.1],...
                                       'callback', {@gaussianMixtureValue, 'a1'}); 
                                   
     S.a2Slider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 6*distance 0.7 0.01],...
               'min',-10,'max',20,...
               'sliderstep',[0.1, 0.1],...
               'callback', {@gaussianMixtureValue, 'a2'});
           
     S.a3Slider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 5*distance 0.7 0.01],...
               'min',-10,'max',20,...
               'sliderstep',[1, 1],...
               'callback', {@gaussianMixtureValue, 'a3'});
           
     S.b1Slider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 3*distance 0.7 0.01],...
               'min',0,'max',20,...
               'sliderstep',[1, 1],...
               'callback', {@gaussianMixtureValue, 'b1'});
           
       S.b2Slider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 2*distance 0.7 0.01],...
               'min',0,'max',20,...
               'sliderstep',[1, 1],...
               'callback', {@gaussianMixtureValue, 'b2'});
           
       S.b3Slider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 1*distance 0.7 0.01],...
               'min',0,'max',20,...
               'sliderstep',[1, 1],...
               'callback', {@gaussianMixtureValue, 'b3'});
           
     guidata(S.fh, S);  

end


function gaussianMixtureValue(hSlider, eventData, param)
S = guidata(hSlider); 
if(param == 'a1')
%     disp("a1 value is ")
%     get(hSlider, 'value')
    S.a1 = get(hSlider, 'value');
elseif(param == 'a2')
%     disp("a2 value is ")
%     get(hSlider, 'value')
    S.a2 = get(hSlider, 'value'); 
elseif(param == 'a3')
    S.a3 = get(hSlider, 'value'); 
elseif(param == 'b1')
    S.b1 = get(hSlider, 'value'); 
elseif(param == 'b2')
    S.b2 = get(hSlider, 'value');
elseif(param == 'b3')
    S.b3 = get(hSlider, 'value'); 
end

x = -5:0.01:5;
z = rectangularPulse(-1, 1, x);
w = 1/exp(-0.5) *  exp(-x.^2/2);

e1=2.9436; e2=15.4697;
f1=30.0000; f2=3.6882;
y2 = e1 * exp(-f1 * x.^2/2)/sqrt(2*pi) + e2* exp(-f2 * x.^2/2)/sqrt(2*pi);

y = S.a1 * exp(-S.b1 * x.^2/2)/sqrt(2*pi) + S.a2 * exp(-S.b2 * x.^2/2)/sqrt(2*pi) ...
    + S.a3 * exp(-S.b3 * x.^2/2)/sqrt(2*pi) ;

plot(x, y, 'g')
hold on
plot(x, z, 'b')
plot(x, w, 'c')
plot(x, y2, 'b', 'LineStyle', ':')
hold off

disp("a value is ")
[S.a1, S.a2, S.a3]
disp("b value is ")
[S.b1, S.b2, S.b3]

guidata(S.fh, S);
end
