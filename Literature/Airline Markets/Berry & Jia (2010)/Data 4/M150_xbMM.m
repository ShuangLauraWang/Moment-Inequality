function [XB alpha p_mc]=M150_xbMM(XMat, theta,modspec,yr)
% This file should accommodate many specifications in both years

fare=XMat(:,1); 
if yr==2006
    switch modspec
        case 20     %demand only
            alpha=theta([1,4]);
            XB=[XMat(:,1:17) * theta([1:3,7:20]), ...
                XMat(:,1:17) * theta([4:6,7:20])];
            p_mc=[];
        case 30     %no LCC
            alpha=theta([1,4]);
            XB=[XMat(:,1:16) * theta([1:3,7:19]), ...
                XMat(:,1:16) * theta([4:6,7:19])];
            p_mc=fare-XMat(:,[17:24,10:16])*theta(20:34);
        case 40     %delay in demand; seat-delay in cost
            alpha=theta([1,4]);
            XB=[XMat(:,1:18) * theta([1:3,7:21]), ...
                XMat(:,1:18) * theta([4:6,7:21])];
            p_mc=fare-XMat(:,[19:28,11:18])*theta(22:39);
        case 42     %seat in cost
            alpha=theta([1,4]);
            XB=[XMat(:,1:17) * theta([1:3,7:20]), ...
                XMat(:,1:17) * theta([4:6,7:20])];
            p_mc=fare-XMat(:,[18:27,10:17])*theta(21:38);
        case 43     %delay in demand; delay in cost
            alpha=theta([1,4]);
            XB=[XMat(:,1:18) * theta([1:3,7:21]), ...
                XMat(:,1:18) * theta([4:6,7:21])];
            p_mc=fare-XMat(:,[19:27,11:18])*theta(22:38);
        case 44     %delay in demand only
            alpha=theta([1,4]);
            XB=[XMat(:,1:18) * theta([1:3,7:21]), ...
                XMat(:,1:18) * theta([4:6,7:21])];
            p_mc=fare-XMat(:,[19:26,11:18])*theta(22:37);
        case 45     %one seat para only
            alpha=theta([1,4]);
            XB=[XMat(:,1:17) * theta([1:3,7:20]), ...
                XMat(:,1:17) * theta([4:6,7:20])];
            p_mc=fare-XMat(:,[18:26,10:17])*theta(21:37);
        case 50     %combine airports; (same as basecase)
            alpha=theta([1,4]);
            XB=[XMat(:,1:17) * theta([1:3,7:20]), ...
                XMat(:,1:17) * theta([4:6,7:20])];
            p_mc=fare-XMat(:,[18:25,10:17])*theta(21:36);
        case 80     %25 airport dummies
            alpha=theta([1,4]);
            XB=[XMat(:,[1:17,26:50]) * theta([1:3,7:20,37:61]), ...
                XMat(:,[1:17,26:50]) * theta([4:6,7:20,37:61])];
            p_mc=fare-XMat(:,[18:25,10:17])*theta(21:36);
        case 90     %Mk>1.5k, 1 set of cost, 1 dist, tour same for both types
            alpha=theta([1,4]);
            XB=[XMat(:,1:16) * theta([1:3,7:19]), ...
                XMat(:,1:16) * theta([4:6,7:19])];
            p_mc=fare-XMat(:,[17:21,9:16])*theta(20:32);
    end;
    
elseif yr==1999
    switch modspec
        case 20 %demand only
            alpha=theta([1,4]);
            XB=[XMat(:,1:18) * theta([1:3,7:21]), ...
                XMat(:,1:18) * theta([4:6,7:21])];
            p_mc=[];
        case 30 %no LCC
            alpha=theta([1,4]);
            XB=[XMat(:,1:18) * theta([1:3,7:21]), ...
                XMat(:,1:18) * theta([4:6,7:21])];
            p_mc=fare-XMat(:,[19:26,10:18])*theta(22:38);
        case 40     %delay in demand; seat-delay in cost
            alpha=theta([1,4]);
            XB=[XMat(:,1:19) * theta([1:3,7:22]), ...
                XMat(:,1:19) * theta([4:6,7:22])];
            p_mc=fare-XMat(:,[20:29,11:19])*theta(23:41);
        case 42     %seat in cost only
            alpha=theta([1,4]);
            XB=[XMat(:,1:18) * theta([1:3,7:21]), ...
                XMat(:,1:18) * theta([4:6,7:21])];
            p_mc=fare-XMat(:,[19:28,10:18])*theta(22:40);
        case 43     %delay in demand; delay in cost
            alpha=theta([1,4]);
            XB=[XMat(:,1:19) * theta([1:3,7:22]), ...
                XMat(:,1:19) * theta([4:6,7:22])];
            p_mc=fare-XMat(:,[20:28,11:19])*theta(23:40);
        case 44     %delay in demand only
            alpha=theta([1,4]);
            XB=[XMat(:,1:19) * theta([1:3,7:22]), ...
                XMat(:,1:19) * theta([4:6,7:22])];
            p_mc=fare-XMat(:,[20:27,11:19])*theta(23:39);
        case 45     %one seat para only
            alpha=theta([1,4]);
            XB=[XMat(:,1:18) * theta([1:3,7:21]), ...
                XMat(:,1:18) * theta([4:6,7:21])];
            p_mc=fare-XMat(:,[19:27,10:18])*theta(22:39);
        case 50     %combine airports; (same as basecase)
            alpha=theta([1,4]);
            XB=[XMat(:,1:18) * theta([1:3,7:21]), ...
                XMat(:,1:18) * theta([4:6,7:21])];
            p_mc=fare-XMat(:,[19:26,10:18])*theta(22:38);
        case 80     %25 airport dummies
            alpha=theta([1,4]);
            XB=[XMat(:,[1:18,27:51]) * theta([1:3,7:21,39:63]), ...
                XMat(:,[1:18,27:51]) * theta([4:6,7:21,39:63])];
            p_mc=fare-XMat(:,[19:26,10:18])*theta(22:38);
        case 90     %Mk>1.5k, 1 set of cost, 1 dist, tour same for both types 
            alpha=theta([1,4]);
            XB=[XMat(:,1:17) * theta([1:3,7:20]), ...
                XMat(:,1:17) * theta([4:6,7:20])];
            p_mc=fare-XMat(:,[18:22,9:17])*theta(21:34);
    end;
end;

