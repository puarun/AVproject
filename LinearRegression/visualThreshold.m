function  pCorrect=visualThreshold(X,Y,T,p,ytest,indext,n)
%   This function will visualize roc curve and effect of picking the threshold
%   Author: Arun.P.U.
    thresh=T(indext);
    Predicted=p;
    Predicted(p>=thresh)=1;
    Predicted(p<thresh)=0;
    [C,order]=confusionmat(ytest,Predicted);
    %pCorrect = (C(1,1)+C(2,2))/(C(1,1)+C(2,2)+C(1,2)+C(2,1));
    pCorrect0 = C(1,1)/(C(1,2)+C(1,1));% TP-FP/TP
    pCorrect1 = C(2,2)/(C(2,1)+C(2,2));
    pCorrect=[pCorrect0,pCorrect1];
    disp(['band ',num2str(n),' Percent correct ',num2str(pCorrect)]); 
    visualize=0;% use to pick visual threshold
    if visualize        
        subplot(2,1,1)
        scatter(p,ytest)
        vline(thresh)
        subplot(2,1,2)
        plot(X,Y),hold on, plot(X(indext),Y(indext),'*'),hold off
        title(['Area under the cuve plot for classifier ',num2str(n)]);
    end
end

