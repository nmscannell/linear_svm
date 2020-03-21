clear all;
load('SVM_data.mat', 'x', 'y');
[m, n] = size(x);
H = zeros(m, m);
for i=1:m
    for j=1:m
        H(i,j) = dot(x(i,:),x(j,:))*y(i)*y(j);
    end
end
f = -(ones(m, 1));
A = y';
b = 0;
l = zeros(m, 1);
alpha = quadprog(H,f,[],[],A,b,l,[]);
mask = alpha>10^(-6);
alpha = alpha.*mask;
indices = (1:m);
svs=indices(alpha>10^(-6))';

w = zeros(m,n);
for i=1:m
    w(i,:) = alpha(i)*y(i)*x(i,:);
end

b = 0;
for i = 1:length(svs)
    b = b+alpha(svs(i))*y(svs(i))*dot(x(svs(i),:),x(27,:));
end
b = 1-b;

figure()
plot(x(1:m/2,1),x(1:m/2,2),'+');
hold on;
plot(x((m/2)+1:m, 1),x((m/2)+1:m,2),'r*');
for i=1:length(svs)
    plot(x(svs(i),1), x(svs(i),2),'ko');
end

x1range = [min(x(:,1))-1, max(x(:,1))+1];
x2range = [min(x(:,2))-1, max(x(:,2))+1];
d = 0.5;
[x1Grid,x2Grid] = meshgrid(x1range(1):d:x1range(2),...
    x2range(1):d:x2range(2));
xGrid = [x1Grid(:),x2Grid(:)];
y_p = sum(w*xGrid')' + b*ones(length(xGrid),1);
contour(x1Grid,x2Grid,reshape(y_p,size(x1Grid)),[0 0],'r','linewidth',2);

sl2 = (x(27,2)-x(26,2))/(x(27,1)-x(26,1));
x0 = linspace(-4,4,100);
y_extra = sl2*x0+(x(27,2)-sl2*x(27,1));
plot(x0, y_extra,'b--');
plot(x0, y_extra+5.35/2, 'k--');
plot(x0, y_extra+5.35,'b--');
ylim([-4, 4]);
xlim([-4, 4]);