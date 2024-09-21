function y2 = yunsuan(a,b,L,t,hfig,h_axes,ci,x3,x4)
%YUNSUAN 此处显示有关此函数的摘要

set(hfig, 'CurrentAxes', h_axes)

%   此处显示详细说明
v=get(t,'value');
Astr = get(a,'String'); %edit_A
A = str2num(Astr);
bstr = get(b,'String'); %edit_B
b = str2num(bstr);
%x=[];%x储存了每次迭代的结果
x11=get(L,'String');
xx=str2num(x11);
x=xx;
% epslion为误差值
epslion=1e-10;
N=100;%最大迭代次数
switch v
    case 1%雅可比迭代法
        L=tril(A,-1); %L为A的单位下三角矩阵
        D=diag(diag(A)); %D为A的对角矩阵
        U=triu(A,1); %U为A的上三角矩阵
        B=-D\(L+U); %B为迭代矩阵
        g=D\b;

        error(1)=1;%用来存放每次误差
        n=1;
        i=1;
        %% 循环
        while i<N
            if error(i)>epslion
                x(:,i+1)=B*x(:,i)+g;

                error(i+1)=norm(A*x(:,i+1)-b);  %可以使用任何范数，此处使用无穷范数
                i=i+1;


                %fprintf('第%d次jacobi迭代结果为：',i-1);
                %disp(x(:,i));
                %fprintf('第%d次迭代结果的误差为：\n',i-1);
                %disp(error(i))

            else
                n=i;
                i=N;
                %fprintf('该方程组的jacobi迭代法的最小迭代次数为：')
                %disp(n-1);
                set(ci,'string',num2str(n-1));
                set(x3,'string',num2str(error));
                set(x4,'string',num2str(x));
            end
        end
        
        j=1;
        Y=[];
        while j<n
            Z(j)=log(error(j+1));
            Y(j)=log(error(j+1));
            j=j+1;

        end
        %plot(X,Y,'b--o','MarkerSize',4*6);
        %h = animatedline(X,Y,'Color','r','LineWidth',1);
        %addpoints(h,X,Y);
        %binscatter(X,Y,20);
        %Z=1:n;
        m=n-1;
        X=2:n;
        plot(X,Y,'.-','MarkerSize',4*6)
        title('雅可比迭代图像','FontSize',12);
        hold on
        p = plot(X(1),Y(1),'o','MarkerFaceColor','red');
        hold off
        axis manual
        for k = 2:length(X)
            p.XData = X(k);
            p.YData = Y(k);
            drawnow
        end
    case 2%getd
        r=b-A*x;				%r为误差
        d=r;					%d为搜索方向
        error(1)=1;%用来存放每次误差
        i=1;
        %% 循环
        while i<N
            if norm(A*x(:,i)-b)>1e-10
                alpha=(r'*d)/(d'*A*d);%步长
                x(:,i+1)=x(:,i)+alpha*d;
                r1=b-A*x(:,i+1);
                error(i+1)=norm(r1);
                bt=-(r1'*A*d)/(d'*A*d);
                d=r1+bt*d;
                r=r1;
                i=i+1;

                %fprintf('第%d次共轭梯度迭代结果为：',i-1);
                %disp(x(:,i));
                %fprintf('第%d次共轭梯度迭代误差为：',i-1);
                %disp(error(i));

            else
                n=i;
                z=x(:,i);
                i=N;
                %fprintf('该方程组的共轭梯度迭代法的最小迭代次数为：')
                %disp(n-1);
                set(ci,'string',num2str(n-1));
                set(x3,'string',num2str(error));
                set(x4,'string',num2str(x));
            end
        end
        j=1;
        while j<n
            Y(j)=log(error(j));
            %Z(j)=log(error(j+1));
            j=j+1;
        end
        m=n-1;
        X=2:n;
        %plot(Y,Z,'-','Linewidth',4*0.5);
        plot(X,Y,'.-','MarkerSize',4*6);
        title('共轭梯度迭代图像','FontSize',12);
        hold on
        p = plot(X(1),Y(1),'o','MarkerFaceColor','red');
        hold off
        axis manual
        for k = 2:length(X)
            p.XData = X(k);
            p.YData = Y(k);
            drawnow
        end
    case 3%SOR
        D=diag(diag(A));%求解对角阵
        LA=tril(A);%抽取下三角矩阵
        U=LA-A;%得到去除对角线元素的上三角阵
        L=D-LA;%得到去除对角线元素的下三角阵
        XXX=[];%存储不同松弛因子的迭代值
        I=[];%记录迭代终止次数
        w=1;%设置松弛因子大小0<w<2
        %% SOR迭代公式
        M=(D-w*L)\((1-w)*D+w*U);
        f=w*((D-w*L)\b);
        MaxEigM = max(abs(eig(M)));%求谱半径
        tu(1)=norm(x);
        i=1;
        while i<N
            if norm(A*x(:,i)-b)>1e-10
                x(:,i+1)=M*x(:,i)+f;
                z=x(:,i+1);
                i=i+1;
                tu(i)=norm(A*x(:,i)-b);
                %fprintf('第%d次SOR迭代结果为：',i-1);%迭代初始值在x第一列
                %disp(x(:,i));
                %fprintf('第%d次SOR迭代误差为：',i-1);
                %disp(norm(A*x(:,i)-b));

            else
                n=i;
                z=x(:,i);
                i=N;
                %fprintf('该方程组的SOR迭代法的最小迭代次数为：')
                %disp(n-1);
                set(ci,'string',num2str(n-1));
                set(x3,'string',num2str(tu));
                set(x4,'string',num2str(x));
            end
        end
        m=n-1;
        X=1:n;
        Y=log(tu(X));
        plot(X,Y,'.-','MarkerSize',4*6);
        title('SOR迭代图像','FontSize',12);
        hold on
        p = plot(X(1),Y(1),'o','MarkerFaceColor','red');
        hold off
        axis manual
        for k = 2:length(X)
            p.XData = X(k);
            p.YData = Y(k);
            drawnow
        end


end

