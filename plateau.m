function plat = plateau(X,opt)
%test = round(abs(randn(100,1)));
X = [0 X.' 0].';
f = find(X == 0);
d = diff(f);
if isempty(f) == 1
    plat = [X (1:length(X)).'];
end

if isempty(f) == 0
if nargin == 1 
    opt = 0;
end

if opt == 0 
    if length(find(d == max(d))) == 1
    plat = [X(f(find(d == max(d)))+1:f(find(d==max(d)))+max(d)-1) (round(f(find(d == max(d)))+1:f(find(d==max(d)))+max(d)-1).')-1];
    else
        %disp('multiple maximum plateaus with same length, outputting full plateau hierarchy')
        opt = 1;
    end
end

if opt == 1
    try
        plat = cell(length(d(d~=1)),1);
        NN = length(d(d~=1));
        for i = unique(sort(d)).'
            if i > 1
                ff = find(d == i);
                N = length(ff);
                while N > 0
                    plat{NN,1} = [X(f(ff(N))+1:f(ff(N))+i-1) (f(ff(N))+1:f(ff(N))+i-1).'];
                    plat{NN,2} = i - 1;
                    NN = NN - 1;
                    N = N - 1;
                end
            end
        end
    catch
        %disp('please input a row vector')
    end
end
end