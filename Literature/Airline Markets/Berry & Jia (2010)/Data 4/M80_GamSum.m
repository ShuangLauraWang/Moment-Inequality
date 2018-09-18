function WMat=M80_GamSum(TMat,gamma,ntype)

if ntype==2
    WMat=TMat(:,1)*gamma(1)+TMat(:,2)*gamma(2);
elseif ntype==3
    WMat=TMat(:,1)*gamma(1)+TMat(:,2)*gamma(2)+TMat(:,3)*gamma(3);
elseif ntype==4
    WMat=TMat(:,1)*gamma(1)+TMat(:,2)*gamma(2)+TMat(:,3)*gamma(3)+TMat(:,4)*gamma(4);
elseif ntype==5
    WMat=TMat(:,1)*gamma(1)+TMat(:,2)*gamma(2)+TMat(:,3)*gamma(3)+TMat(:,4)*gamma(4)+...
        TMat(:,5)*gamma(5);
end;
