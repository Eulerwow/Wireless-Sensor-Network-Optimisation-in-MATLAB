%% plot the location of sensors, one of intitial point of the relay, and the opimised location of the relay

scatter(X(:,1),X(:,2),'filled','square','blue');
hold on
scatter(s0(:,1),s0(:,2),'filled','diamond','black');

scatter(MLSmins(:,1),MLSmins(:,2),'filled','diamond','green');
scatter(L2smins(:,1),L2smins(:,2),'filled','diamond','red');