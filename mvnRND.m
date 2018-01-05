function para = mvnRND(muv,mun,sigv,sign,p)
x = 0:0.01:10;
y = 1*((abs(x)./sigv.^2).*exp(-((x-muv).^2)./(2*sigv.^2)));
Y = y;
Y = round(Y*1000);
Y(1) = 1;
vel = zeros(sum(Y),1);
for i = 1:length(Y)-1
    vel(sum(Y(1:i)):sum(Y(1:i+1))) = x(i);
end
n = 1;
Z = mvnrnd([0 0], [1 p; p 1], n);
U = normcdf(Z);
para = [quantile(vel,U(:,1)) abs(norminv(U(:,2),mun,sign))];