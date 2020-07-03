function dydt = chen(t, y)
    dydt = zeros(3, 1);
    dydt(1) = 35*(y(2) -y(1));
    dydt(2) = -7 * y(1) - y(1) * y(3) + 28 * y(2);
    dydt(3) = y(1) * y(2) - y(3);
end