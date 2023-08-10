function [rmse,Sk] = get_Skill(pred,obs)
% This function

yhat = pred;
y = obs;

rmse = sqrt(mean((y - yhat).^2));  % Root Mean Squared Error

Sk = 1-sum((y-yhat).^2)./sum((y-mean(y)).^2);
