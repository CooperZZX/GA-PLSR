%不同求导方法对比
clear
close all
restoredefaultpath;
load('LSCC-data.mat');
order=2;
framelen=5;
[~,g] = sgolay(order,framelen);
spec_diff=diff(LSCC_Spec_Sele(1,:))./0.005;
d_spec0 = conv(LSCC_Spec_Sele(1,:),  g(:,1), 'same'); 
d_spec1 = conv(LSCC_Spec_Sele(1,:),  1/(-0.005) * g(:,2)); 
plot(LSCC_Spec_Sele(1,:));
hold on
plot(d_spec0);
y = sgolayfilt(LSCC_Spec_Sele(1,:),order,framelen);
plot(y);
figure();
plot(spec_diff);
hold on
plot(d_spec1);
plot(diff(d_spec0)./0.005);