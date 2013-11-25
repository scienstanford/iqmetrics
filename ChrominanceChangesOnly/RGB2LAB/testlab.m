RGB1 = double(imread('hats.jpg'));
% figure;
% imshow(RGB1);
% convert to XYZ
XYZ1 = rgb2xyz(RGB1);
% convert to Lab
lab1 = xyz2lab(XYZ1);

% change   L
 L1 = lab1(:,:,1);
 l1 = 0.2*(L1 - mean(L1(:))) + mean(L1(:));
 a1 = lab1(:,:,2);
 b1 = lab1(:,:,3);
% figure; hist(L1(:));title('L1');
% figure; hist(l1(:));title('l1')
% figure; plot(L1(:),l1(:));
% xlabel('L1');
% ylabel('l1');



 XYZ2= lab2xyz(lab1);
 RGB2 = xyz2rgb(XYZ2);
 
%  figure;
%  imshow(RGB1/255);
%  figure;
%  imshow(RGB2/255);
 
 XYZ3 = rgb2xyz(RGB2*255);
lab2 = xyz2lab(XYZ3);

 L2 = lab2(:,:,1);
 a2 = lab1(:,:,2);
 b2 = lab1(:,:,3);


% figure;
% subplot(2,2,1),plot(lab1(:,:,1),lab2(:,:,1),'*');title('L');xlabel('L1');ylabel('L2');
% subplot(2,2,2),plot(lab1(:,:,2),lab2(:,:,2),'*');title('a');xlabel('a1');ylabel('a2');
% subplot(2,2,3),plot(lab1(:,:,3),lab2(:,:,3),'*');title('b');xlabel('b1');ylabel('b2');

% LAB with gray(luminance);

RGB3 = rgb2gray(RGB1/255);
% figure;
% subplot(2,2,1),plot(lab1(:,:,1),RGB3,'*');title('L');xlabel('L1');ylabel('gray');
% subplot(2,2,2),plot(lab1(:,:,2),RGB3,'*');title('a');xlabel('a1');ylabel('gray');
% subplot(2,2,3),plot(lab1(:,:,3),RGB3,'*');title('b');xlabel('b1');ylabel('gray');


% exchange a and b;

lab3 = lab1;
lab3(:,:,2)= b1;
lab3(:,:,3)= a1;

 XYZ4= lab2xyz(lab3);
 RGB4 = xyz2rgb(XYZ2);
 
 XYZ5 = rgb2xyz(RGB4*255);
lab5 = xyz2lab(XYZ4);

 L3 = lab5(:,:,1);
 a3 = lab5(:,:,2);
 b3 = lab5(:,:,3);
 
 XYZ6= lab2xyz(lab5);
 RGB6 = xyz2rgb(XYZ6);
 
%  figure;
%  imshow(RGB6);
 
 %  figure the exchange a and b
% figure;
% subplot(2,2,1),plot(lab1(:,:,1),lab5(:,:,1),'*');title('L');xlabel('L1');ylabel('L5');
% subplot(2,2,2),plot(lab1(:,:,2),lab5(:,:,2),'*');title('a');xlabel('a1');ylabel('a5');
% subplot(2,2,3),plot(lab1(:,:,3),lab5(:,:,3),'*');title('b');xlabel('b1');ylabel('b5');

%%%%

figure;

       subplot(3,3,1);imshow(a1);title('Original Image'); 
       H = fspecial('motion',9,0);
       MotionBlura1 = imfilter(a1,H,'replicate');
       subplot(3,3,2);imshow(MotionBlura1);title('Motion Blurred Image1');
       H = fspecial('motion',20,0);
       MotionBlura2 = imfilter(a1,H,'replicate');
       subplot(3,3,3);imshow( MotionBlura2);title('Motion Blurred Image2');
       H = fspecial('motion',30,0);
       MotionBlura3 = imfilter(a1,H,'replicate');
       subplot(3,3,3);imshow(MotionBlura3);title('Motion Blurred Image3');
       H = fspecial('motion',40,0);
       MotionBlura4 = imfilter(a1,H,'replicate');
       subplot(3,3,4);imshow(MotionBlura4);title('Motion Blurred Image3');
       H = fspecial('motion',50,0);
       MotionBlura5 = imfilter(a1,H,'replicate');
       subplot(3,3,5);imshow(MotionBlura5);title('Motion Blurred Image3');

% figure;
%        subplot(3,3,1);imshow(a1);title('Original Image'); 
%        H = fspecial('motion',20,45);
%        MotionBlura = imfilter(a1,H,'replicate');
%        subplot(3,3,2);imshow(MotionBlura);title('Motion Blurred Image');
%        H = fspecial('disk',10);
%        blurreda = imfilter(a1,H,'replicate');
%        subplot(3,3,3);imshow(blurreda);title('Blurred Image');
%        H = fspecial('unsharp');
%        sharpeneda = imfilter(a1,H,'replicate');
%        subplot(3,3,4);imshow(sharpeneda);title('Sharpened Image');
%        Spnoisea = imnoise(a1,'salt & pepper', 0.02);
%        subplot(3,3,5);imshow(Spnoisea);title('salt & pepper');
%        gaussiana = imnoise(a1,'gaussian', 0.02);
%        subplot(3,3,6);imshow(gaussiana);title('gaussian');
%        specklea = imnoise(a1,'speckle', 0.02);
%        subplot(3,3,7);imshow(specklea);title('speckle');
  
%   produce different image.

            lab6(:,:,1)= L3;
            lab6(:,:,2) = MotionBlura1;
            lab6(:,:,3) = b3;
            
            lab7(:,:,1)= L3;
            lab7(:,:,2) = MotionBlura2;
            lab7(:,:,3) = b3;
            
            lab8(:,:,1)= L3;
            lab8(:,:,2) = MotionBlura3;
            lab8(:,:,3) = b3;
            
            lab9(:,:,1)= L3;
            lab9(:,:,2) = MotionBlura4;
            lab9(:,:,3) = b3;
            
            
             XYZ4a= lab2xyz(lab6);
            RGB4a = xyz2rgb(XYZ4a);
            
             figure;
             imshow(RGB4a/255);
             
           XYZ4a = rgb2xyz(RGB4a);

           lab4a = xyz2lab(XYZ4a); 
           
            L4a = lab4a(:,:,1);
           
            RGB4a = rgb2gray(RGB4a/255);
            
            
              XYZ7a= lab2xyz(lab7);
            RGB7a = xyz2rgb(XYZ7a);
            
             figure;
             imshow(RGB7a/255);
             
             
           XYZ7a = rgb2xyz(RGB7a);

           lab7a = xyz2lab(XYZ7a); 
           
           L7a = lab7a(:,:,1);
           
            RGB7a = rgb2gray(RGB7a/255)
            
            
                     XYZ8a= lab2xyz(lab8);
            RGB8a = xyz2rgb(XYZ8a);
            
             figure;
             imshow(RGB8a/255);
             
           XYZ8a = rgb2xyz(RGB8a);

           lab8a = xyz2lab(XYZ8a); 
           
           L8a = lab8a(:,:,1);
           
            RGB8a = rgb2gray(RGB8a/255)
            
            
            
                     XYZ9a= lab2xyz(lab9);
            RGB9a = xyz2rgb(XYZ9a);
            
             figure;
             imshow(RGB9a/255);
             
           XYZ9a = rgb2xyz(RGB9a);

           lab9a = xyz2lab(XYZ9a); 
           
           L9a = lab9a(:,:,1);
           
            RGB9a = rgb2gray(RGB9a/255);
            
            
               
            
            figure;
            subplot(3,2,1),plot(lab1(:,:,1),RGB3,'*');title('before blur');xlabel('L1');ylabel('gray');
            subplot(3,2,2),plot(L4a,RGB4a,'*');title('After blur ');xlabel('L1');ylabel('gray');
            subplot(3,2,3),plot(L7a,RGB7a,'*');title('After blur ');xlabel('L1');ylabel('gray');
            subplot(3,2,4),plot(L8a,RGB8a,'*');title('After blur ');xlabel('L1');ylabel('gray');
            subplot(3,2,5),plot(L9a,RGB9a,'*');title('After blur ');xlabel('L1');ylabel('gray');
            %subplot(3,2,6),plot(L4a,RGB4a,'*');title('After blur ');xlabel('L1');ylabel('gray');
       
       
  figure;
       subplot(3,2,1),plot(a1,MotionBlura,'*');title('MotionBlur');xlabel('a1');ylabel('MotionBlur');
       subplot(3,2,2),plot(a1,blurreda,'*');title('blurreda');xlabel('a1');ylabel('blurred');
       subplot(3,2,3),plot(a1,sharpeneda,'*');title('Sharpened');xlabel('a1');ylabel('Sharpened');
       subplot(3,2,4),plot(a1,Spnoisea,'*');title('salt & pepper');xlabel('a1');ylabel('salt & pepper');    
        subplot(3,2,5),plot(a1,gaussiana,'*');title('gaussian');xlabel('a1');ylabel('gaussian');
       subplot(3,2,6),plot(a1,specklea,'*');title('speckle');xlabel('a1');ylabel('speckle'); 
        
        %%%%%
      figure;title('a channel'); 
      subplot(3,3,1);imshow(b1);title('Original Image'); 
       H = fspecial('motion',20,45);
       MotionBlurb = imfilter(b1,H,'replicate');
       subplot(3,3,2);imshow(MotionBlurb);title('Motion Blurred Image');
       H = fspecial('disk',10);
       blurredb = imfilter(b1,H,'replicate');
       subplot(3,3,3);imshow(blurredb);title('Blurred Image');
       H = fspecial('unsharp');
       sharpenedb = imfilter(b1,H,'replicate');
       subplot(3,3,4);imshow(sharpenedb);title('Sharpened Image');
       Spnoiseb = imnoise(b1,'salt & pepper', 0.02);
       subplot(3,3,5);imshow(Spnoisea);title('salt & pepper');
       gaussianb = imnoise(b1,'gaussian', 0.02);
       subplot(3,3,6);imshow(gaussianb);title('gaussian');
       speckleb = imnoise(b1,'speckle', 0.02);
       subplot(3,3,7);imshow(speckleb);title('speckle');
       
       figure;
       subplot(3,2,1),plot(b1,MotionBlurb,'*');title('MotionBlur');xlabel('b1');ylabel('MotionBlur');
       subplot(3,2,2),plot(b1,blurredb,'*');title('blurreda');xlabel('b1');ylabel('blurred');
        subplot(3,2,3),plot(b1,sharpenedb,'*');title('Sharpened');xlabel('b1');ylabel('Sharpened');
        subplot(3,2,4),plot(b1,Spnoiseb,'*');title('salt & pepper');xlabel('b1');ylabel('salt & pepper');    
        subplot(3,2,5),plot(b1,gaussianb,'*');title('gaussian');xlabel('b1');ylabel('gaussian');
       subplot(3,2,6),plot(b1,speckleb,'*');title('speckle');xlabel('b1');ylabel('speckle'); 


