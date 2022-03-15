% simple nls
clear all
close all
t = [1; 2; 4; 5; 8];
y = [3.2939; 4.2699; 7.1749; 9.3008; 20.259];


x_1 = 20;
x_2 = 0.7;

theta = [x_1; x_2];



J = [exp(theta(2)*t(1)),t(1)*theta(1)*exp(theta(2)*t(1));
    exp(theta(2)*t(2)),t(2)*theta(1)*exp(theta(2)*t(2));
    exp(theta(2)*t(3)),t(3)*theta(1)*exp(theta(2)*t(3));
    exp(theta(2)*t(4)),t(4)*theta(1)*exp(theta(2)*t(4));
    exp(theta(2)*t(5)),t(5)*theta(1)*exp(theta(2)*t(5));];


for i=1:height(t)
    error(i,1) = y(i)-theta(1)*exp(theta(2)*t(i));
end


while(norm(error)>9e-5)
figure(1)
clf
plot(t,y)
hold on
plot(t, theta(1)*exp(theta(2)*t))



 norm(error)

 delta_y = inv(J'*J)*J'*error;
 theta = theta+delta_y;

J = [exp(theta(2)*t(1)),t(1)*theta(1)*exp(theta(2)*t(1));
    exp(theta(2)*t(2)),t(2)*theta(1)*exp(theta(2)*t(2));
    exp(theta(2)*t(3)),t(3)*theta(1)*exp(theta(2)*t(3));
    exp(theta(2)*t(4)),t(4)*theta(1)*exp(theta(2)*t(4));
    exp(theta(2)*t(5)),t(5)*theta(1)*exp(theta(2)*t(5));];





for i=1:height(t)
    error(i,1) = y(i)-theta(1)*exp(theta(2)*t(i));
end

end