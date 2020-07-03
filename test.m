
img = imread('lena256.bmp');
Q_R = img(:,:,1);
Q_G = img(:,:,2);
Q_B = img(:,:,3);

NN = 5000
[M,N]=size(Q_R);
x1=ceil(rand(1,NN)*(M-1));      %����5000��1~M-1�����������Ϊ��
y1=ceil(rand(1,NN)*(N-1));      %����5000��1~N-1�����������Ϊ��

%ˮƽ
XX_R_SP=zeros(1,NN);YY_R_SP=zeros(1,NN);  %Ԥ�����ڴ�
XX_G_SP=zeros(1,NN);YY_G_SP=zeros(1,NN);
XX_B_SP=zeros(1,NN);YY_B_SP=zeros(1,NN);
%��ֱ
XX_R_CZ=zeros(1,NN);YY_R_CZ=zeros(1,NN);  %Ԥ�����ڴ�
XX_G_CZ=zeros(1,NN);YY_G_CZ=zeros(1,NN);
XX_B_CZ=zeros(1,NN);YY_B_CZ=zeros(1,NN);
%�Խ���
XX_R_DJX=zeros(1,NN);YY_R_DJX=zeros(1,NN);  %Ԥ�����ڴ�
XX_G_DJX=zeros(1,NN);YY_G_DJX=zeros(1,NN);
XX_B_DJX=zeros(1,NN);YY_B_DJX=zeros(1,NN);
for i=1:NN
    %ˮƽ
    XX_R_SP(i)=Q_R(x1(i),y1(i));
    YY_R_SP(i)=Q_R(x1(i)+1,y1(i));
    XX_G_SP(i)=Q_G(x1(i),y1(i));
    YY_G_SP(i)=Q_G(x1(i)+1,y1(i));
    XX_B_SP(i)=Q_B(x1(i),y1(i));
    YY_B_SP(i)=Q_B(x1(i)+1,y1(i));
    %��ֱ
    XX_R_CZ(i)=Q_R(x1(i),y1(i));
    YY_R_CZ(i)=Q_R(x1(i),y1(i)+1);
    XX_G_CZ(i)=Q_G(x1(i),y1(i));
    YY_G_CZ(i)=Q_G(x1(i),y1(i)+1);
    XX_B_CZ(i)=Q_B(x1(i),y1(i));
    YY_B_CZ(i)=Q_B(x1(i),y1(i)+1);
    %�Խ���
    XX_R_DJX(i)=Q_R(x1(i),y1(i));
    YY_R_DJX(i)=Q_R(x1(i)+1,y1(i)+1);
    XX_G_DJX(i)=Q_G(x1(i),y1(i));
    YY_G_DJX(i)=Q_G(x1(i)+1,y1(i)+1);
    XX_B_DJX(i)=Q_B(x1(i),y1(i));
    YY_B_DJX(i)=Q_B(x1(i)+1,y1(i)+1);
end
%ˮƽ
figure;scatter(XX_R_SP,YY_R_SP,18,'filled');title('ԭͼRͨ��ˮƽ����Ԫ�������');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
figure;scatter(XX_G_SP,YY_G_SP,18,'filled');title('ԭͼGͨ��ˮƽ����Ԫ�������');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
figure;scatter(XX_B_SP,YY_B_SP,18,'filled');title('ԭͼBͨ��ˮƽ����Ԫ�������');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
%��ֱ
figure;scatter(XX_R_CZ,YY_R_CZ,18,'filled');title('ԭͼRͨ����ֱ����Ԫ�������');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
figure;scatter(XX_G_CZ,YY_G_CZ,18,'filled');title('ԭͼGͨ����ֱ����Ԫ�������');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
figure;scatter(XX_B_CZ,YY_B_CZ,18,'filled');title('ԭͼBͨ����ֱ����Ԫ�������');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
%�Խ���
figure;scatter(XX_R_DJX,YY_R_DJX,18,'filled');title('ԭͼRͨ���Խ�������Ԫ�������');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
figure;scatter(XX_G_DJX,YY_G_DJX,18,'filled');title('ԭͼGͨ���Խ�������Ԫ�������');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
figure;scatter(XX_B_DJX,YY_B_DJX,18,'filled');title('ԭͼBͨ���Խ�������Ԫ�������');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
%Rͨ��
Q_R=double(Q_R);
EX2_R=0;EY2_SP_R=0;DX2_R=0;DY2_SP_R=0;COVXY2_SP_R=0;    %ˮƽ
EY2_CZ_R=0;DY2_CZ_R=0;COVXY2_CZ_R=0;    %��ֱ
EY2_DJX_R=0;DY2_DJX_R=0;COVXY2_DJX_R=0;   %�Խ���
%Gͨ��
Q_G=double(Q_G);
EX2_G=0;EY2_SP_G=0;DX2_G=0;DY2_SP_G=0;COVXY2_SP_G=0;    %ˮƽ
EY2_CZ_G=0;DY2_CZ_G=0;COVXY2_CZ_G=0;    %��ֱ
EY2_DJX_G=0;DY2_DJX_G=0;COVXY2_DJX_G=0;   %�Խ���
%Bͨ��
Q_B=double(Q_B);
EX2_B=0;EY2_SP_B=0;DX2_B=0;DY2_SP_B=0;COVXY2_SP_B=0;    %ˮƽ
EY2_CZ_B=0;DY2_CZ_B=0;COVXY2_CZ_B=0;    %��ֱ
EY2_DJX_B=0;DY2_DJX_B=0;COVXY2_DJX_B=0;   %�Խ���
for i=1:NN
    %��һ�����ص��E��ˮƽ����ֱ���Խ���ʱ����ó��ĵ�һ�����ص��E��ͬ��ͳһ��EX2��ʾ
    EX2_R=EX2_R+Q_R(x1(i),y1(i));
    EX2_G=EX2_G+Q_G(x1(i),y1(i));
    EX2_B=EX2_B+Q_B(x1(i),y1(i));
    %�ڶ������ص��E��ˮƽ����ֱ���Խ��ߵ�E�ֱ��ӦEY2_SP��EY2_CZ��EY2_DJX
    %Rͨ��
    EY2_SP_R=EY2_SP_R+Q_R(x1(i),y1(i)+1);
    EY2_CZ_R=EY2_CZ_R+Q_R(x1(i)+1,y1(i));
    EY2_DJX_R=EY2_DJX_R+Q_R(x1(i)+1,y1(i)+1);
    %Gͨ��
    EY2_SP_G=EY2_SP_G+Q_G(x1(i),y1(i)+1);
    EY2_CZ_G=EY2_CZ_G+Q_G(x1(i)+1,y1(i));
    EY2_DJX_G=EY2_DJX_G+Q_G(x1(i)+1,y1(i)+1);
    %Bͨ��
    EY2_SP_B=EY2_SP_B+Q_B(x1(i),y1(i)+1);
    EY2_CZ_B=EY2_CZ_B+Q_B(x1(i)+1,y1(i));
    EY2_DJX_B=EY2_DJX_B+Q_B(x1(i)+1,y1(i)+1);
end
%ͳһ��ѭ����������ص����1000���ɼ����������
%Rͨ��
EX2_R=EX2_R/NN;
EY2_SP_R=EY2_SP_R/NN;
EY2_CZ_R=EY2_CZ_R/NN;
EY2_DJX_R=EY2_DJX_R/NN;
%Gͨ��
EX2_G=EX2_G/NN;
EY2_SP_G=EY2_SP_G/NN;
EY2_CZ_G=EY2_CZ_G/NN;
EY2_DJX_G=EY2_DJX_G/NN;
%Bͨ��
EX2_B=EX2_B/NN;
EY2_SP_B=EY2_SP_B/NN;
EY2_CZ_B=EY2_CZ_B/NN;
EY2_DJX_B=EY2_DJX_B/NN;

for i=1:NN
    %��һ�����ص��D��ˮƽ����ֱ���Խ���ʱ����ó���һ�����ص��D��ͬ��ͳһ��DX2��ʾ
    DX2_R=DX2_R+(Q_R(x1(i),y1(i))-EX2_R)^2;
    DX2_G=DX2_G+(Q_G(x1(i),y1(i))-EX2_G)^2;
    DX2_B=DX2_B+(Q_B(x1(i),y1(i))-EX2_B)^2;
    %�ڶ������ص��D��ˮƽ����ֱ���Խ��ߵ�D�ֱ��ӦDY2_SP��DY2_CZ��DY2_DJX
    %Rͨ��
    DY2_SP_R=DY2_SP_R+(Q_R(x1(i),y1(i)+1)-EY2_SP_R)^2;
    DY2_CZ_R=DY2_CZ_R+(Q_R(x1(i)+1,y1(i))-EY2_CZ_R)^2;
    DY2_DJX_R=DY2_DJX_R+(Q_R(x1(i)+1,y1(i)+1)-EY2_DJX_R)^2;
    %Gͨ��
    DY2_SP_G=DY2_SP_G+(Q_G(x1(i),y1(i)+1)-EY2_SP_G)^2;
    DY2_CZ_G=DY2_CZ_G+(Q_G(x1(i)+1,y1(i))-EY2_CZ_G)^2;
    DY2_DJX_G=DY2_DJX_G+(Q_G(x1(i)+1,y1(i)+1)-EY2_DJX_G)^2;
    %Bͨ��
    DY2_SP_B=DY2_SP_B+(Q_B(x1(i),y1(i)+1)-EY2_SP_B)^2;
    DY2_CZ_B=DY2_CZ_B+(Q_B(x1(i)+1,y1(i))-EY2_CZ_B)^2;
    DY2_DJX_B=DY2_DJX_B+(Q_B(x1(i)+1,y1(i)+1)-EY2_DJX_B)^2;
    %�����������ص���غ����ļ��㣬ˮƽ����ֱ���Խ���
    %Rͨ��
    COVXY2_SP_R=COVXY2_SP_R+(Q_R(x1(i),y1(i))-EX2_R)*(Q_R(x1(i),y1(i)+1)-EY2_SP_R);
    COVXY2_CZ_R=COVXY2_CZ_R+(Q_R(x1(i),y1(i))-EX2_R)*(Q_R(x1(i)+1,y1(i))-EY2_CZ_R);
    COVXY2_DJX_R=COVXY2_DJX_R+(Q_R(x1(i),y1(i))-EX2_R)*(Q_R(x1(i)+1,y1(i)+1)-EY2_DJX_R);
    %Gͨ��
    COVXY2_SP_G=COVXY2_SP_G+(Q_G(x1(i),y1(i))-EX2_G)*(Q_G(x1(i),y1(i)+1)-EY2_SP_G);
    COVXY2_CZ_G=COVXY2_CZ_G+(Q_G(x1(i),y1(i))-EX2_G)*(Q_G(x1(i)+1,y1(i))-EY2_CZ_G);
    COVXY2_DJX_G=COVXY2_DJX_G+(Q_G(x1(i),y1(i))-EX2_G)*(Q_G(x1(i)+1,y1(i)+1)-EY2_DJX_G);
    %Bͨ��
    COVXY2_SP_B=COVXY2_SP_B+(Q_B(x1(i),y1(i))-EX2_B)*(Q_B(x1(i),y1(i)+1)-EY2_SP_B);
    COVXY2_CZ_B=COVXY2_CZ_B+(Q_B(x1(i),y1(i))-EX2_B)*(Q_B(x1(i)+1,y1(i))-EY2_CZ_B);
    COVXY2_DJX_B=COVXY2_DJX_B+(Q_B(x1(i),y1(i))-EX2_B)*(Q_B(x1(i)+1,y1(i)+1)-EY2_DJX_B);
end
%ͳһ��ѭ����������ص����1000���ɼ����������
%Rͨ��
DX2_R=DX2_R/NN;
DY2_SP_R=DY2_SP_R/NN;
DY2_CZ_R=DY2_CZ_R/NN;
DY2_DJX_R=DY2_DJX_R/NN;
COVXY2_SP_R=COVXY2_SP_R/NN;
COVXY2_CZ_R=COVXY2_CZ_R/NN;
COVXY2_DJX_R=COVXY2_DJX_R/NN;
%Gͨ��
DX2_G=DX2_G/NN;
DY2_SP_G=DY2_SP_G/NN;
DY2_CZ_G=DY2_CZ_G/NN;
DY2_DJX_G=DY2_DJX_G/NN;
COVXY2_SP_G=COVXY2_SP_G/NN;
COVXY2_CZ_G=COVXY2_CZ_G/NN;
COVXY2_DJX_G=COVXY2_DJX_G/NN;
%Bͨ��
DX2_B=DX2_B/NN;
DY2_SP_B=DY2_SP_B/NN;
DY2_CZ_B=DY2_CZ_B/NN;
DY2_DJX_B=DY2_DJX_B/NN;
COVXY2_SP_B=COVXY2_SP_B/NN;
COVXY2_CZ_B=COVXY2_CZ_B/NN;
COVXY2_DJX_B=COVXY2_DJX_B/NN;
%ˮƽ����ֱ���Խ��ߵ������
%Rͨ��
RXY2_SP_R=COVXY2_SP_R/sqrt(DX2_R*DY2_SP_R);
RXY2_CZ_R=COVXY2_CZ_R/sqrt(DX2_R*DY2_CZ_R);
RXY2_DJX_R=COVXY2_DJX_R/sqrt(DX2_R*DY2_DJX_R);
%Gͨ��
RXY2_SP_G=COVXY2_SP_G/sqrt(DX2_G*DY2_SP_G);
RXY2_CZ_G=COVXY2_CZ_G/sqrt(DX2_G*DY2_CZ_G);
RXY2_DJX_G=COVXY2_DJX_G/sqrt(DX2_G*DY2_DJX_G);
%Bͨ��
RXY2_SP_B=COVXY2_SP_B/sqrt(DX2_B*DY2_SP_B);
RXY2_CZ_B=COVXY2_CZ_B/sqrt(DX2_B*DY2_CZ_B);
RXY2_DJX_B=COVXY2_DJX_B/sqrt(DX2_B*DY2_DJX_B);

%% ���������Ϣ
disp('Rͨ������ԣ�');
disp(['��ԭͼƬRͨ������ԣ�','  ˮƽ�����=',num2str(RXY2_SP_R),'  ��ֱ�����=',num2str(RXY2_CZ_R),'  �Խ��������=',num2str(RXY2_DJX_R)]);
disp('Gͨ������ԣ�');
disp(['��ԭͼƬGͨ������ԣ�','  ˮƽ�����=',num2str(RXY2_SP_G),'  ��ֱ�����=',num2str(RXY2_CZ_G),'  �Խ��������=',num2str(RXY2_DJX_G)]);
disp('Bͨ������ԣ�');
disp(['��ԭͼƬBͨ������ԣ�','  ˮƽ�����=',num2str(RXY2_SP_B),'  ��ֱ�����=',num2str(RXY2_CZ_B),'  �Խ��������=',num2str(RXY2_DJX_B)]);
