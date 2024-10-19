function P = objectiveFcn(C)
    [v_f, v_z, t]=get_v(C);%调用函数
    P =sum(abs(v_f-v_z).^2)*C*(1/20);
end

%将浮子和振子的速度储存在[v_f, v_z, t]=get_v(C)这个函数里面

function [v_f, v_z, t]=get_v(C)
    tspan = [180,200];%180s到200s状态稳定
    y0 = [0 0 0 0];
    [t, x] = ode45(@(t, x) ode1(t, x, C), tspan, y0);
    v_f = x(:, 2);
    v_z = x(:, 4);
end

%用龙格库塔算法计算浮子和振子的位移和速度
function dx=ode1(t, x, C)
    m = 4866;
    m_z = 2433;
    A = 1165.992;
    B = 167.8395;
    K = 80000;
    e = 1025;
    g = 9.8;
    S = pi;
    f = 4890;
    w = 2.2143;
    dx = zeros(4, 1);
    dx(1) = x(2);
    dx(2) = (-B*x(2)+C*(x(4)-x(2))+K*(x(3)-x(1))-(e*g*S)*x(1)+f*cos(w*t))/(m+A);
    dx(3) = x(4);
    dx(4) = (-C*(x(4)-x(2))- K*(x(3)-x(1)))/m_z;
end