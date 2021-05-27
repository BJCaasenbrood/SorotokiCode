%clr;
load('fem18May');

%% make shape libary
shp = Shapes([0,1,0,0,0,0],Y,'NModal',8);

shp = shp.fitPOD();
shp = shp.rebuild('NModal',1);

% for ii = 1:1:size(shp.posData,2)
%     N = shp.posData{ii};
%     z(ii,:) = N(1,:);
% end
% p = shp.string(325);
% plot(p(:,1),p(:,3)); axis equal; hold on;
% plot(N(:,1),N(:,2)); 
% 
% Q = linspace(0,337,25);
% Ba = shp.get('Ba');
% XI{1} = zeros(80,1);
% 
% for ii = 1:1:numel(Q)
% [P,~,S] = shp.string(Q(ii));
% 
% ee = zeros(1,length(S));
% for jj = 1:length(S)
%     ee(:,jj) = shp.Phi(S(jj))*Q(ii);
% end
% 
% XI{1} = vappend(XI{1},ee(1,:).',2);
% end

% X = linspace(0,1,80);
% % S = shp.xiData{1};
% figure(103); cla; subplot(2,2,2);
% %Q = linspace(0,337,25);
% XII = XI{1};
% for ii = 1:1:26
% 
% %     P = shp.string(Q(ii));
% %     P = P*1e3;
%     plot(S/0.12,XII(:,ii),'Linewidth',1.5,...
%         'Color',col(1,ii/(1.05*26)));
%      hold on;
% %     axis equal;
% %     z2(ii,:) = [P(end,1),P(end,3)];
%     %Y{ii} = 1e-3*[Nx,Ny];
% end

% plot(z(:,1)*1e3,z(:,2)*1e3,'Linewidth',1.5,...
%         'Color',col(1)); hold on;
% plot(z2(:,1),z2(:,2),'Linewidth',1.5,...
%         'Color',col(2)); hold on;    
    
% XI = shp.xiData{1};
% C = XI*XI.';
% 
% [U,V,D] = svd(C);
% % 
% figure(101)
% for ii = 1:26
%     plot(XI(:,ii));
%     hold on;
% end
% 
% figure(102)
% subplot(2,2,1);
% semilogy(diag(V))
% 
% subplot(2,2,2);
% for ii = 1:9
%     plot(U(:,ii));
%     hold on;
% end
