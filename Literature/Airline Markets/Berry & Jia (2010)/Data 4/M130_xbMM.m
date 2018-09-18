function [XB alpha p_mc]=M130_xbMM(XMat, theta,modspec,yr)
% This file should accommodate many specifications in both years

fare=XMat(:,1); 
if yr==2006
    switch modspec
        case 22     %2 types, 3 coef (include a type-specific const), 1 dept,
            alpha=theta([1,4]);
            XB=[XMat(:,1:17) * theta([1:3,7:20]), ...
                XMat(:,1:17) * theta([4:6,7:20])];
            p_mc=fare-XMat(:,[3,18:21,10:17])*theta(21:33);
        case 68     %allow diff cost coef for short-medium/long haul flights
            alpha=theta([1,4]);
            XB=[XMat(:,1:17) * theta([1:3,7:20]), ...
                XMat(:,1:17) * theta([4:6,7:20])];
            p_mc=fare-XMat(:,[18:25,10:17])*theta(21:36);
        case 69     %2 sets of cost coef; tour=0 for business type
            alpha=theta([1,5]);
            XB=[XMat(:,1:17) * theta([1:4,8:20]), ...
                XMat(:,[1:2,4:17]) * theta([5:7,8:20])];
            p_mc=fare-XMat(:,[18:25,10:17])*theta(21:36);
        case 70     %2 sets of cost coef; tour coef type specific
            alpha=theta([1,5]);
            XB=[XMat(:,1:17) * theta([1:4,9:21]), ...
                XMat(:,1:17) * theta(5:21)];
            p_mc=fare-XMat(:,[18:25,10:17])*theta(22:37);
        case 71     %demand only; 2 sets of cost coef; tour=0 for business type
            alpha=theta([1,5]);
            XB=[XMat(:,1:17) * theta([1:4,8:20]), ...
                XMat(:,[1:2,4:17]) * theta([5:7,8:20])];
            p_mc=[];
        case 72     %2 sets of cost coef; tour=0 for business type; no B6
            alpha=theta([1,5]);
            XB=[XMat(:,1:16) * theta([1:4,8:19]), ...
                XMat(:,[1:2,4:16]) * theta([5:7,8:19])];
            p_mc=fare-XMat(:,[17:24,10:16])*theta(20:34);
        case 73     %delay in demand; seat/delay in supply
            alpha=theta([1,5]);
            XB=[XMat(:,1:18) * theta([1:4,8:21]), ...
                XMat(:,[1:2,4:18]) * theta([5:7,8:21])];
            p_mc=fare-XMat(:,[19:28,11:18])*theta(22:39);
        case 74     % 25 airport dummies
            alpha=theta([1,5]);
            XB=[XMat(:,[1:17,26:50]) * theta([1:4,8:20,37:61]), ...
                XMat(:,[1:2,4:17,26:50]) * theta([5:7,8:20,37:61])];
            p_mc=fare-XMat(:,[18:25,10:17])*theta(21:36);
        case 75     %Mk>3k,tour=0 for busi
            alpha=theta([1,5]);
            XB=[XMat(:,1:17) * theta([1:4,8:20]), ...
                XMat(:,[1:2,4:17]) * theta([5:7,8:20])];
            p_mc=fare-XMat(:,[18:22,10:17])*theta(21:33);
        case 76     %Mk>3k, no B6/SW Entry, tour=0 for busi
            alpha=theta([1,5]);
            XB=[XMat(:,1:16) * theta([1:4,8:19]), ...
                XMat(:,[1:2,4:16]) * theta([5:7,8:19])];
            p_mc=fare-XMat(:,[17:21,10:16])*theta(20:31);
        case 77     %Mk>3k,2 tour coef
            alpha=theta([1,5]);
            XB=[XMat(:,1:17) * theta([1:4,9:21]), ...
                XMat(:,1:17) * theta([5:8,9:21])];
            p_mc=fare-XMat(:,[18:22,10:17])*theta(22:34);
        case 78     %Mk>3k,1 tour coef, 1 dist
            alpha=theta([1,4]);
            XB=[XMat(:,1:16) * theta([1:3,7:19]), ...
                XMat(:,1:16) * theta([4:6,7:19])];
            p_mc=fare-XMat(:,[17:21,9:16])*theta(20:32);
    end;
    
elseif yr==1999
    switch modspec
        case 22 %2 types, 3 coef (include a type-specific const), 1 dept
            alpha=theta([1,4]);
            XB=[XMat(:,1:18) * theta([1:3,7:21]), ...
                XMat(:,1:18) * theta([4:6,7:21])];
            p_mc=fare-XMat(:,[3,19:22,10:18])*theta(22:35);
        case 68 %cost coef differ for short_medium/long haul routes
            alpha=theta([1,4]);
            XB=[XMat(:,1:18) * theta([1:3,7:21]), ...
                XMat(:,1:18) * theta([4:6,7:21])];
            p_mc=fare-XMat(:,[19:26,10:18])*theta(22:38);
        case 69 %cost coef differ for short_medium/long haul routes
            alpha=theta([1,5]);
            XB=[XMat(:,1:18) * theta([1:4,8:21]), ...
                XMat(:,[1:2,4:18]) * theta(5:21)];
            p_mc=fare-XMat(:,[19:26,10:18])*theta(22:38);
        case 70 %2 set of cost coef, tour coef type specific
            alpha=theta([1,5]);
            XB=[XMat(:,1:18) * theta([1:4,9:22]), ...
                XMat(:,1:18) * theta(5:22)];
            p_mc=fare-XMat(:,[19:26,10:18])*theta(23:39);
        case 71 %demand only; 2 set of cost coef, tour coef type specific
            alpha=theta([1,5]);
            XB=[XMat(:,1:18) * theta([1:4,8:21]), ...
                XMat(:,[1:2,4:18]) * theta(5:21)];
            p_mc=[];
        case 73 %delay in demand; delay/seat in supply
            alpha=theta([1,5]);
            XB=[XMat(:,1:19) * theta([1:4,8:22]), ...
                XMat(:,[1:2,4:19]) * theta(5:22)];
            p_mc=fare-XMat(:,[20:29,11:19])*theta(23:41);
        case 74 %25 airport dummies
            alpha=theta([1,5]);
            XB=[XMat(:,[1:18,27:51]) * theta([1:4,8:21,39:63]), ...
                XMat(:,[1:2,4:18,27:51]) * theta([5:21,39:63])];
            p_mc=fare-XMat(:,[19:26,10:18])*theta(22:38);
        case 75 % Mk longer than 3k
            alpha=theta([1,5]);
            XB=[XMat(:,1:18) * theta([1:4,8:21]), ...
                XMat(:,[1:2,4:18]) * theta(5:21)];
            p_mc=fare-XMat(:,[19:23,10:18])*theta(22:35);
        case 77     %Mk>3k,2 tour coef
            alpha=theta([1,5]);
            XB=[XMat(:,1:18) * theta([1:4,9:22]), ...
                XMat(:,1:18) * theta([5:8,9:22])];
            p_mc=fare-XMat(:,[19:23,10:18])*theta(23:36);
        case 78     %Mk>3k,1 tour coef, 1 dist
            alpha=theta([1,4]);
            XB=[XMat(:,1:17) * theta([1:3,7:20]), ...
                XMat(:,1:17) * theta([4:6,7:20])];
            p_mc=fare-XMat(:,[18:22,9:17])*theta(21:34);
    end;
end;

