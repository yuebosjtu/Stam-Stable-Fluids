clear;
% writerObj=VideoWriter('F:/paperResults/databox/boxdorp1.mp4','MPEG-4');  
 writerObj=VideoWriter('pond_TurJet_WNEO.mp4','MPEG-4');  
%  writerObj=VideoWriter('blowbox/blowboxs.mp4','MPEG-4');  
% Set frame rate
writerObj.FrameRate = 30;
% Open video writer object and write frames sequentially
open(writerObj)
re = zeros(256,512*3,3);
for i = 1:1:146        % Some number of frames
     % Read frame
      
%      frame = sprintf('blowbox/smoke/ball_data%04d.ppm', i);
%     frame = sprintf('F:/paperResults/databox/smoke_box/ball_data%04d.ppm', i);
%     frame1 = sprintf('img/dambreak%05d.png', i);
%     frame2 = sprintf('ppm_ve_r/im%05d.ppm', i);
%     frame1 = sprintf('ppm_phi_Y/static%05d.ppm', i);
%     frame2 = sprintf('ppm_Pressure_Y/static%05d.ppm', i);
%    frame3 = sprintf('ppm_ve_Y/static%05d.ppm', i);

 frame1 = sprintf('ppm_ve/im%05d.ppm', i);
    frame2 = sprintf('ppm_rho/im%05d.ppm', i);
   frame3 = sprintf('ppm_t/im%05d.ppm', i);   
%    frame4 = sprintf('ppm_Pressure/im%05d.ppm', i);
      input1 = imread(frame1);
      input2 = imread(frame2);
       input3 = imread(frame3); 
%             input4 = imread(frame4); 
        re(:,1:512,:)= input1;
         re(:,512+1:512*2,:)= input2;
          re(:,512*2+1:512*3,:)= input3;
%         re(:,256*3+1:256*4,:)= input4; 
     re = uint8(re);
     writeVideo(writerObj, re);  
   
%         frame = sprintf('renderBasketHigh/flag%05d.ppm', i);     
%       input1 = imread(frame);
%        writeVideo(writerObj, input1);
end
%
% Close the video writer objectsdad
close(writerObj);
