clear all
%构建界面
set(gcf,'DefaultUicontrolUnits', 'Normal',...
    'Name','不同迭代法求解方程组Ax=b', 'NumberTitle','off',...
    'MenuBar','none') ;
%text部分
uicontrol(gcf,'style','text', 'position',[0.05,0.9,0.07,0.046],...
    'string','输入A', 'fontsize',12,...
    'BackgroundColor',[0.75 0.75 0.75]); 
uicontrol(gcf,'style','text', 'position',[0.05,0.65,0.07,0.046],...
    'string','输入b', 'fontsize',12,...
    'BackgroundColor',[0.75 0.75 0.75]); 
uicontrol(gcf,'style','text', 'position',[0.05,0.25,0.08,0.046],...
    'string','初始值', 'fontsize',12,...
    'BackgroundColor',[0.75 0.75 0.75]); 
uicontrol(gcf,'style','text', 'position',[0.05,0.125,0.1,0.046],...
    'string','初始误差', 'fontsize',12,...
    'BackgroundColor',[0.75 0.75 0.75]); 
uicontrol(gcf,'style','text', 'position',[0.55,0.85,0.1,0.05],...
    'string','迭代方式', 'fontsize',12,...
    'BackgroundColor',[0.75 0.75 0.75]); 
uicontrol(gcf,'style','text', 'position',[0.55,0.35,0.1,0.05],...
    'string','迭代次数', 'fontsize',12,...
    'BackgroundColor',[0.75 0.75 0.75]); 
uicontrol(gcf,'style','text', 'position',[0.55,0.25,0.03,0.05],...
    'string','第', 'fontsize',12,...
    'BackgroundColor',[0.75 0.75 0.75]);
uicontrol(gcf,'style','text', 'position',[0.69,0.25,0.12,0.05],...
    'string','次迭代误差', 'fontsize',12,...
    'BackgroundColor',[0.75 0.75 0.75]);
uicontrol(gcf,'style','text', 'position',[0.55,0.15,0.06,0.05],...
    'string','误差:', 'fontsize',12,...
    'BackgroundColor',[0.75 0.75 0.75]);
uicontrol(gcf,'style','text', 'position',[0.74,0.15,0.05,0.05],...
    'string','x值:', 'fontsize',12,...
    'BackgroundColor',[0.75 0.75 0.75]);

%输入
a = uicontrol(gcf, 'Style','edit',...
    'FontSize',12,'Pos',[0.14,0.78,0.335,0.18],...
  'max', 3, 'min',0, 'HorizontalAlignment','Left');
b = uicontrol(gcf, 'Style','edit',...
    'FontSize',12,'Pos',[0.14,0.53,0.335,0.18],...
   'max', 3, 'min',0, 'HorizontalAlignment','Left');
ic = uicontrol(gcf, 'Style','edit',...
    'FontSize',12,'Pos',[0.59,0.25,0.09,0.05],...
   'max', 3, 'min',0, 'HorizontalAlignment','Left');

%选择迭代法
t = uicontrol(gcf,'Style','popupmenu', ...
    'String', ...
    '雅可比迭代法|共轭梯度法|SOR迭代法', ...
    'Position',[0.65,0.7,0.25,0.2]);



%结果
L = uicontrol(gcf, 'Style','Edit',...
    'FontSize',12,'Pos',[0.18, 0.225,0.275,0.1],...
    'max', 3, 'min',0,'HorizontalAlignment','Left');
U = uicontrol(gcf, 'Style','Edit',...
    'FontSize',12,'Pos',[0.18, 0.1 ,0.275,0.1],...
    'max', 3, 'min',0,'HorizontalAlignment','Left');
x3 = uicontrol(gcf, 'Style','Edit',...
    'FontSize',12,'Pos',[423 398 100 23],...
    'max', 3, 'min',0,'HorizontalAlignment','Left');
x4 = uicontrol(gcf, 'Style','Edit',...
    'FontSize',12,'Pos',[423 390 100 23],...
    'max', 3, 'min',0,'HorizontalAlignment','Left');
ci = uicontrol(gcf, 'Style','Edit',...
    'FontSize',12,'Pos',[0.7,0.35,0.1,0.05],...
    'HorizontalAlignment','Left');
wu = uicontrol(gcf, 'Style','Edit',...
    'FontSize',12,'Pos',[0.63,0.15,0.09,0.05],...
    'max', 3, 'min',0, ...
    'HorizontalAlignment','Left');
xz = uicontrol(gcf, 'Style','Edit',...
    'FontSize',12,'Pos',[0.81,0.15,0.09,0.05],...
    'max', 3, 'min',0, ...
    'HorizontalAlignment','Left');

%运算
push_A = uicontrol(gcf,'style','push', 'string','随机生成初始值',...
    'pos',[0.2,0.4,0.22,0.05],'fontsize',12, ...
    'call',...
    'y1 = suiji(a,b,L,U);'); 

%push_B = uicontrol(gcf,'style','push','String','确认', ...程滢晓3019210107
   % 'Position',[0.7,0.75,0.1,0.05],'fontsize',12, ...
   % 'call', ...
    %'y2 = yunsuan(a,b,t,L);2,');程滢晓3019210107
push_draw=uicontrol(gcf, 'Style','Push',...
    'String','确认','value',0,'FontSize',12,......
    'Pos',[0.7,0.75,0.1,0.05], ...
    'Call', 'yunsuan(a,b,L,t,gcf,h_axes,ci,x3,x4)');

push_C = uicontrol(gcf,'style','push', 'string','确认',...
    'pos',[0.83,0.25,0.06,0.05],'fontsize',12, ...
    'call',...
    'y3 = diicijieguo(wu,xz,ic,x3,x4);'); 

%图片

h_axes = axes('Position',[0.6,0.5,0.3,0.2]);

title('迭代图像','FontSize',12);