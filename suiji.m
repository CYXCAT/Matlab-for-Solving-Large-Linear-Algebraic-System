function y1 = suiji(a,b,L,U)
%SUIJI 此处显示有关此函数的摘要
%   此处显示详细说明
Astr = get(a,'String'); %edit_A
A = str2num(Astr);
bstr = get(b,'String'); %edit_B
b = str2num(bstr);
x=[];%x储存了每次迭代的结果
x(:,1)=5.*randn(size(b));%给出任意的迭代初始值，占了第一列
d=b-A*x;
r=norm(d);
set(L,'string',num2str(x)); 
set(U,'string',num2str(r));
end

