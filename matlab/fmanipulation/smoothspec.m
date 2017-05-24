function x = smoothspec(x, n)

    w = hann(n);
    w = w/sum(w);
    for i=1:size(x,2)
        xx = x(:,i);
        xx = conv2([ones(n,1)*xx(1); xx(:); ones(n,1)*xx(end)], w, 'same');
        xx = xx(n+1:end-n);
        x(:,i) = xx;
    end
end
