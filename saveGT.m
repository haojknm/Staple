% % % gtdata = [];
% % % for i = 1:223
% % %     num = 2524+(i-1)*2+100000;
% % %     str = num2str(num);
% % %     I = imread( [ 'E:\Project\tracking\����\staple-master\sequences\person1\img_' str(2:end) '.bmp' ] );
% % %     imshow(I)
% % %     hold on
% % %     x=[];
% % %     y=[];
% % %     n=0;
% % %     disp('������������ȡ��Ҫ�ĵ�');
% % %     disp('��������Ҽ���ȡ���һ����');
% % %     but=1;
% % % 
% % %     for j=1:2
% % %     [xi,yi,but]=ginput(1);
% % %     plot(xi,yi,'bo')
% % %     x(j,1)=xi;
% % %     y(j,1)=yi;
% % %     end
% % %     hold off
% % %     close
% % %     
% % %     gtdata(i,3) = x(2) - x(1) + 1;
% % %     gtdata(i,4) = y(2) - y(1) + 1;
% % %     gtdata(i,1) = x(1) -1 + gtdata(i,3)/2;
% % %     gtdata(i,2) = y(1) -1 + gtdata(i,4)/2;
% % %     
% % %     str = num2str(num+1);
% % %     I2 = imread( [ 'G:\dataset\�������ݿ�\Dataset 03 OSU Color and Thermal Database\frames\6a\img_' str(2:end) '.bmp' ] );
% % %     str = num2str(i+100000);
% % %     imwrite( I2, [ 'E:\Project\tracking\����\staple-master\sequences\person1\imgs\person1' str(2:end) '.bmp']);
% % % end
% % % save person1 gtdata




clear
gtdata = [];
for i = 1:541
    num = 508+(i-1)*2+100000;
    str = num2str(num);
    I = imread( [ 'E:\Project\tracking\����\staple-master\sequences\person2\bk\img_' str(2:end) '.bmp' ] );
    imshow(I)
    set(gcf,'Position',get(0,'ScreenSize'))
    hold on
    x=[];
    y=[];
    n=0;
    disp('������������ȡ��Ҫ�ĵ�');
    disp('��������Ҽ���ȡ���һ����');
    but=1;

    for j=1:2
    [xi,yi,but]=ginput(1);
    plot(xi,yi,'bo')
    x(j,1)=xi;
    y(j,1)=yi;
    end
    hold off
    close
    
    gtdata(i,3) = x(2) - x(1) + 1;
    gtdata(i,4) = y(2) - y(1) + 1;
    gtdata(i,1) = x(1) -1 + gtdata(i,3)/2;
    gtdata(i,2) = y(1) -1 + gtdata(i,4)/2;
    
    str = num2str(num+1);
    I2 = imread( [ 'G:\dataset\�������ݿ�\Dataset 03 OSU Color and Thermal Database\frames\5a\img_' str(2:end) '.bmp' ] );
    str = num2str(i+100000);
    imwrite( I2, [ 'E:\Project\tracking\����\staple-master\sequences\person2\imgs\person2' str(2:end) '.bmp']);
end
save person2 gtdata