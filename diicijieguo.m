function y3 = diicijieguo(wu,xz,ic,x3,x4)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
icstr = get(ic,'String'); 
ic = str2num(icstr);
i=ic+1;
x3str = get(x3,'String'); 
x3 = str2num(x3str);
x4str = get(x4,'String'); 
x = str2num(x4str);
xn = x(:,i);
a=x3(i);
set(wu,'string',num2str(a));
set(xz,'string',num2str(xn));
end

