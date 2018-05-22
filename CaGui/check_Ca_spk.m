myfigure('CheckSpkCa');clf;
ii = 62;
plot(Ff(ii,:),'b');

spkind = find(Spk(ii,:));
line(repmat(spkind,2,1),(max(ylim)-min(ylim))/10*[full(Spk(ii,spkind)); zeros(1,length(spkind))]+min(Ff(ii,:)),'Color',[1 0 0]);