function [idx, cm2] = LgsTs(l, u, x0, N)
    cm = zeros(1, l*2+N);
    cm(1) = x0;
    for i = 2:(l*2+N)
        cm(i) = u * cm(i-1) * (1-cm(i-1));
    end
    cm1 = cm(N+1:l+N);
    cm2 = cm(l+N+1:end);
    [~, idx] = sort(cm1);
end
