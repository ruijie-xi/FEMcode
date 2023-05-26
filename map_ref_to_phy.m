function pts = map_ref_to_phy(vertices,ref_pts)

x1=vertices(1,1);
y1=vertices(2,1);
x2=vertices(1,2);
y2=vertices(2,2);
x3=vertices(1,3);
y3=vertices(2,3);

pts = ref_pts;
pts(1,:)=x1+(x2-x1)*ref_pts(1,:)+(x3-x1)*ref_pts(2,:);
pts(2,:)=y1+(y2-y1)*ref_pts(1,:)+(y3-y1)*ref_pts(2,:);
end