function finaldec = D2A_Conversion_Quant(bitstream ,n)

scalefactor = 2^(n-1);

finaldec = zeros(round(numel(bitstream)/n),1);
f=1;fg=1;
qnt_val=0;
while f <= numel(bitstream) - (n-1)
    for i=0:n-1
    qnt_val = qnt_val + (2^(n-1-i))*bitstream(f+i);
    end
    decimal_value = qnt_val ;
    finaldec(fg)=decimal_value;
    f=f+n;
    fg=fg+1;
    qnt_val=0;
end

finaldec = finaldec/(0.5*scalefactor);
mxf = max(finaldec);% mnf = min(finaldec);
finaldec = finaldec/mxf;
finaldec = finaldec - ones(numel(finaldec),1);

end