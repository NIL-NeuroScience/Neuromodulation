function downsampled = f_downsample(sig,bin)

if bin == 1
    downsampled = sig;
else
    [H,W,~] = size(sig);

    H_trim = floor(H / bin) * bin;
    W_trim = floor(W / bin) * bin;
    sig = sig(1:H_trim,1:W_trim,:);

    sig = movmean(sig,bin,1,EndPoints='discard');
    sig = movmean(sig,bin,2,EndPoints='discard');

    downsampled = sig(1:bin:end,1:bin:end,:);
end
