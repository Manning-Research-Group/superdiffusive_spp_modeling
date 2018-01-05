function [thresh] = otsu (c,plotToggle)
% changing from [s,u,thresh] to [thresh]
% Shellhammer-Pavlidis algorithm 
% "Selective sampling" + Otsu's algorithm

% Sample at peaks: Assume c does not start or end with a peak.
[T,temp] = size(c);
s_index = 1;
for c_index = 2:T-1
  if c(c_index-1) < c(c_index) & c(c_index+1) <= c(c_index)  % local max
    s(s_index) = c(c_index);
	s_index = s_index + 1;
  elseif c(c_index-1) >= c(c_index) & c(c_index+1) > c(c_index)  % local max
    s(s_index) = c(c_index);
	s_index = s_index + 1;
  end;
end;

% Otsu's Algorithm
best_var = inf;
thresh = 0.01;
sorted_s = sort(abs(s));
[temp,S] = size(sorted_s);
for i = 1:S-1
  new_var = cov(sorted_s(1:i)) + cov(sorted_s(i+1:S));
  if new_var < best_var
    best_var = new_var;
	thresh = (sorted_s(i)+sorted_s(i+1)) / 2;
  end;
end;
    
% Threshold s: prune out all values < thresh
u_index = 1;
for s_index = 1:S
  if abs(s(s_index)) >= thresh
    u(u_index) = s(s_index);
	u_index = u_index + 1;
  end;
end;

% Plots
if plotToggle == 1
subplot(4,1,1);
plot(c);
title('Original derivative signal');

subplot(4,1,2);
stem(s);
ylabel('Sampling');

subplot(4,1,3);
hist(abs(s));
hold on;
a = axis;
x = [thresh,thresh];
y = [a(3),a(4)];
plot(x,y,'r');
ylabel('Histogram');
hold off;

subplot(4,1,4);
stem(s);
[N,temp] = size(s);
hold on;
a = axis;
x = [a(1),a(2)];
y = [thresh,thresh];
plot(x,y,'r');
plot(x,-y,'r');
ylabel(strcat('Thresh = ',num2str(thresh)));
hold off;
end
% subplot(5,1,5);
% plot(u);
% ylabel('Final');
