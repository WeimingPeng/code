a = func1(pi*2*5,2,5)
function func1 = func1(v_sep,z,l)
    % 计算液位高度的函数，本函数将多次用到
    % z为积分起始 l为传递面长度
    syms h
    a = v_sep/l;
    r = 2;
    y = sqrt(2*r*h - h^2); % 注意此处r应于par.r相同
    c = 2*int(y,z,h);
    eq_sep1 = c-a;
    yy = solve(eq_sep1,h);
    h = real(yy);
    func1 = h;
end