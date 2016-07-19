function xhamp = run_hampel(x)
xhamp = zeros(size(x,1), size(x,2));
    for i=1:size(x,2)
        Gene = x(:,i);
        Genehamp = hampel(Gene, size(x,1), 3);
        xhamp(:,i) = zscore(Genehamp);
    end