function out = Encoder(bits)

    out=zeros(round(numel(bits)/2),1);
    m=1;
    k=1;
    while k <= numel(bits)-2
        num = (bits(k+1)) + 2*(bits(k));
        out(m)=num;
        m=m+1;
        k=k+2;
    end

end