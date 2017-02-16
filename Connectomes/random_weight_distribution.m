function [suffled, values]=random_weight_distribution(Adj)

[str] = strengths_und(Adj); 
[deg] = degrees_und(Adj);
values=cell(length(deg),1);
    for i=1:length(deg)
            r=rand(1,deg(i));
            S=sum(r);
            y=r/S;
            value=y*str(i);
            values{i}=value; 
    end
    
flag=1;
X=find(Adj);
suffled=zeros(length(X),1);
        for i=1:length(deg)
    m=length(values{i});
    suffled(flag:flag+m-1) = values{i};
    flag = flag+m;
        end
end