for qq = 1:10:n
%     figure;
    
    Ez(20,:,qq)
    imagesc(Ez(:,:,qq));
    xlabel(qq);
    colormap jet
    colorbar
    caxis([-0.1,0.1]);
    pause(0.3);
%     close all
%     
%     
% %     plot([1:size(Ey,1)],Ey(:,qq),'r');
% %     xlabel(['t=' num2str(delta_t*(qq-1)*1e9) 'ns']);
% %     axis([0 size(Ey,1) 1.1*min(Ey,[],'ALL') 1.1*max(Ey,[],'ALL')]);
end