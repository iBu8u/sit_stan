function [model] = RevLearn_RL_beta_alt2(params)

%% prepare data and initilize parameters

load('../_data/data3_129.mat')
load('../_outputs/param_alt2_MAT.mat')


lr    = params(:,1:2);
evidW = params(:,3:4);
bta   = params(:,6:10);

ns = size(data3,3);
nt = size(data3,1);

c1 = squeeze(data3(:,3,:))';
c2 = squeeze(data3(:,10,:))';
r  = squeeze(data3(:,14,:))';

otherC1   = zeros(ns,nt,4);
wOthers   = zeros(ns,nt,4);
cfsC2     = zeros(ns,nt,4);
cfoC2     = zeros(ns,nt,4);

prob_sC2  = zeros(ns,nt,4);
prob_oC2  = zeros(ns,nt,4);
pOth_C2   = zeros(ns,nt,4);

prob_sC2_ori  = zeros(ns,nt,4);
prob_oC2_ori  = zeros(ns,nt,4);
pOth_C2_ori   = zeros(ns,nt,4);

wProb_sC2 = zeros(ns,nt,4);
wProb_sC2 = zeros(ns,nt,4);
with    = zeros(ns,nt);
against = zeros(ns,nt);

for s = 1:ns
    otherC1(s,:,:) = data3(:,6:9,s);
    wOthers(s,:,:) = data3(:,51:54,s);
    cfsC2(s,:,:)   = data3(:,61:64,s);
    cfoC2(s,:,:)   = data3(:,65:68,s);
    prob_sC2_ori(s,:,:) = data3(:,73:76,s);
    prob_oC2_ori(s,:,:) = data3(:,77:80,s);
    
    for t = 1:nt
        my1    = data3(t,3,s);
        other1 = data3(t,6:9,s);
        with(s,t)    = length(find(other1==my1));
        against(s,t) = length(find(other1~=my1));
    end
end

%% model

for s = 1:ns
    myV  = zeros(nt+1,2);
    othV = zeros(nt+1,2);
    prob = zeros(nt+1,2);
    pe   = zeros(nt,1);
    penc = zeros(nt,1);
    
    for t = 1:nt
        for o = 1:4
            prob_sC2(s,t,o)  = betacdf(0.5, evidW(s,1)*cfoC2(s,t,o)+1, evidW(s,2)*cfsC2(s,t,o)+1 );
            prob_oC2(s,t,o)  = 1 - prob_sC2(s,t,o);
            wProb_sC2(s,t,o) = wOthers(s,t,o) * prob_sC2(o);
            wProb_oC2(s,t,o) = wOthers(s,t,o) * prob_oC2(o);
        end
        
        if c2(s,t) == 1
            pOth_C2(s,t,:) = prob_sC2(s,t,:);
            pOth_C2_ori(s,t,:) = prob_sC2_ori(s,t,:);
        elseif c2(s,t) == 2
            pOth_C2(s,t,:) = prob_oC2(s,t,:);
            pOth_C2_ori(s,t,:) = prob_oC2_ori(s,t,:);
        end
    end
end

pOth_C2 = pOth_C2(:,1:99,:);
pOth_C2_ori = pOth_C2_ori(:,1:99,:);
otherC1 = otherC1(:,2:100,:);
otherC1(otherC1==2) = 0;

cm = zeros(ns,4); % correlation matrix

for s = 1:ns
    
    a  = squeeze(pOth_C2(s,:,:))>0.5;
    b  = squeeze(pOth_C2_ori(s,:,:))>0.5;
    c  = squeeze(otherC1(s,:,:));
    
    for o= 1:4
        % acc(s,o) = sum(sum([a(:,o) b(:,o)],2)==2) / length(a);
        acc1(s,o) = sum(a(:,o)==c(:,o)) / length(a);
        acc2(s,o) = sum(b(:,o)==c(:,o)) / length(b);
        
    end
    
%     An=bsxfun(@minus,a,mean(a,1));
%     Bn=bsxfun(@minus,b,mean(b,1));
%     An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));
%     Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));
%     cm(s,:)=sum(An.*Bn,1);
    
end



%% output
model.prob_sC2 = prob_sC2;





return;
