function vflat=LinMapsqueeze(sflat,par)
s=reshape(sflat,par.Ndim,[]);
v=LinMap(s,par);
vflat=v(:);
end