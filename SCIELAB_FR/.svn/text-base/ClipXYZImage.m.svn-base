function resImage=ClipXYZImage(xyzImage, whitePt)
% function dstImage=ClipXYZImage(xyzImage, whitePt)
% Purpose: This function applies clipping to the incoming image "xyzImage".
% The resultant image will be clipped into within the range [0,
% white_point]. 

[M, N, L]=size(xyzImage);

for ii=1:L
  t = xyzImage(:,:,ii);
  t(t<0) = 0; t(t>whitePt(ii)) = whitePt(ii);
  resImage(:,:,ii) = t;
end
