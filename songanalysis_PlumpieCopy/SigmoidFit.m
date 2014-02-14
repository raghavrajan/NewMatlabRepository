function error = SigmoidFit(init, x, y)

error = sum((init(3) + init(4)*(1./(1+exp((init(1)-x)*init(2)))) - y).^2);