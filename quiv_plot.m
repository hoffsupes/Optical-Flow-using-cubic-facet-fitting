
function quiv_plot(U,V,n)       

sacl = n/size(U,2);         %% scaling factor
unew = sacl*imresize(U,sacl);  %% new U and V
vnew = sacl*imresize(V,sacl);
quiver(unew(end:-1:1,:),-vnew(end:-1:1,:),2); %% quiver because it's coordinate system opposite to matlab's
axis('tight');
end