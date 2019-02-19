function varargout = lapulassp(varargin)
% LAPULASSP MATLAB code for lapulassp.fig
%      LAPULASSP, by itself, creates a new LAPULASSP or raises the existing
%      singleton*.
%
%      H = LAPULASSP returns the handle to a new LAPULASSP or the handle to
%      the existing singleton*.
%
%      LAPULASSP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LAPULASSP.M with the given input arguments.
%
%      LAPULASSP('Property','Value',...) creates a new LAPULASSP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lapulassp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lapulassp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lapulassp

% Last Modified by GUIDE v2.5 13-Dec-2018 11:45:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lapulassp_OpeningFcn, ...
                   'gui_OutputFcn',  @lapulassp_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%RP金字塔

function Y = RP_fuse(M1, M2, zt, ap, mp)
%Y = fuse_rat(M1, M2, zt, ap, mp) image fusion with ratio pyramid
%
%    M1 - input image A
%    M2 - input image B
%    zt - maximum decomposition level
%    ap - coefficient selection highpass (see selc.m) 
%    mp - coefficient selection base image (see selb.m) 
%
%    Y  - fused image   

%    (Oliver Rockinger 16.08.99)

% check inputs 
[z1,s1] = size(M1);
[z2,s2] = size(M2);
if (z1 ~= z2) || (s1 ~= s2)
  error('Input images are not of same size');
end

% define filter 
w  = [1 4 6 4 1] / 16;

% define eps
eps = 1e-6;

% cells for selected images
E = cell(1,zt);

% loop over decomposition depth -> analysis
for i1 = 1:zt 
  % calculate and store actual image size 
  [z,s]  = size(M1); 
  zl(i1) = z; sl(i1)  = s;
  
  % check if image expansion necessary 
  if (floor(z/2) ~= z/2), ew(1) = 1; else, ew(1) = 0; end
  if (floor(s/2) ~= s/2), ew(2) = 1; else, ew(2) = 0; end

  % perform expansion if necessary
  if (any(ew))
  	M1 = adb(M1,ew);
  	M2 = adb(M2,ew);
  end
  
  % perform filtering 
  G1 = conv2(conv2(es2(M1,2), w, 'valid'),w', 'valid');
  G2 = conv2(conv2(es2(M2,2), w, 'valid'),w', 'valid');
 
  % decimate, undecimate and interpolate 
  M1T = conv2(conv2(es2(undec2(dec2(G1)), 2), 2*w, 'valid'),2*w', 'valid');
  M2T = conv2(conv2(es2(undec2(dec2(G2)), 2), 2*w, 'valid'),2*w', 'valid');
 
  % select coefficients and store them
  E(i1) = {selc(M1./(M1T+eps), M2./(M2T+eps), ap)};
  
  % decimate 
  M1 = dec2(G1);
  M2 = dec2(G2);
end

% select base coefficients of last decompostion stage
M1 = selb(M1,M2,mp);

% loop over decomposition depth -> synthesis
for i1 = zt:-1:1
  % undecimate and interpolate 
  M1T = conv2(conv2(es2(undec2(M1), 2), 2*w, 'valid'), 2*w', 'valid');
  % add coefficients
  M1  = (M1T+eps) .* E{i1};
  % select valid image region 
  M1 	= M1(1:zl(i1),1:sl(i1));
end

% copy image
Y = M1;

%小波变换wt

function y= wtfusion(x1,x2,N,wname)

%函数功能：
%     函数x=wtfusion(x1,x2,N,wname)将两幅原图像x1,x2进行基于小波变换的图像融合，得到融合后的图像y
%     近似分量采用加权平均的融合规则，各细节分量采用基于区域特性量测的融合规则
%输入参数：
%     x1----输入原图像1
%     x2----输入原图像2
%     N----小波分解的层数
%     wname----小波基函数
%输出参数：
%     y----原图像融合后得到的图像
%-----------------------------------------------------------------%

x1=double(x1);                   %将uint8的图像数据类型转换成double型进行数据处理
x2=double(x2);

 %将原图像x1,x2分别进行N层小波分解，wname为小波基函数，
 %C为各层分解系数,S为各层分解系数长度,也就是大小. 
 %C的结构:c=[A(N)|H(N)|V(N)|D(N)|H(N-1)|V(N-1)|D(N-1)|H(N-2)|V(N-2)|D(N-2)|...|H(1)|V(1)|D(1)]
 %A(N)代表第N层低频系数,H(N)|V(N)|D(N)代表第N层高频系数,分别是水平,垂直,对角高频
 %
 %S是一个(N+2)行、2列的结构是储存各层分解系数长度的,即第一行是A(N)的长度（其实是A(N)的原矩阵的行数和列数）,
 %第二行是H(N)|V(N)|D(N)|的长度,第三行是H(N-1)|V(N-1)|D(N-1)的长度,
 %倒数第二行是H(1)|V(1)|D(1)长度,最后一行是X的长度(大小)

[C1,S1]=wavedec2(x1,N,wname);    % 函数功能：实现图像的小波分解。wname是小波基的名称。N是分解的层数.
[C2,S2]=wavedec2(x2,N,wname);    

A1=appcoef2(C1,S1,wname,N);            %函数功能：离散小波变换低频部分系数提取
A2=appcoef2(C2,S2,wname,N);             %提取出小波分解的近似分量。低频只有一个，就是A(N)。其余的都是高频。
A=0.5*A1+0.5*A2;                       %低频分量的融合规则采用加权平均的方法

%figure
%imshow(uint8(A))

%仿照matlab中近似分量和细节分量的存储方式，把融合后的近似分量和细节分量转成行向量，然后存入向量C中
%这样做是为了方便重构原图像

a=reshape(A,1,S1(1,1)*S1(1,2));        %将A转换成行向量.
C=a;

for i=N:-1:1                           %循环从第N层到第1层    
    [H1,V1,D1]=detcoef2('all',C1,S1,i);       %小波变换高频部分提取
    [H2,V2,D2]=detcoef2('all',C2,S2,i);        %提取出小波分解的各层细节分量
    H=f(H1,H2);                           
    V=f(V1,V2);    %高频分量的融合规则：采用基于区域特性量测的融合规则。可以选择不同的融合规则。
    D=f(D1,D2);
    h=reshape(H,1,S1(N+2-i,1)*S1(N+2-i,2));%分别将融合后的细节分量转成行向量，并存入行向量C中
    v=reshape(V,1,S1(N+2-i,1)*S1(N+2-i,2));
    d=reshape(D,1,S1(N+2-i,1)*S1(N+2-i,2));
    C=[C,h,v,d];
end

S=S1;
y=waverec2(C,S,wname);      %重构原图像
%figure(1);imshow(uint8(y));title('基于小波变换的融合图像')


%小波变换f子函数
function y = f(x1 , x2)

%函数功能：
%       y=f(x1,x2)将两幅原图像x1和x2基于区域特性量测的融合规则进行融合,得到融合后的图像y
%       首先计算两幅图像的匹配度，若匹配度大于阈值，说明两幅图像对应局部能量较接近，
%       因此采用加权平均的融合方法；若匹配度小于阈值，说明两幅图像对应局部能量相差较大，
%       因此选取局部区域能量较大的小波系数作为融合图像的小波系数
%输入参数：
%      x1----输入原图像1
%      x2----输入原图像2
%输出参数：
%      y----融合后的图像
%------------------------------------------------------------%

w=1/16*[1 2 1;2 4 2;1 2 1];   %权系数
E1=conv2(x1.^2,w,'same');     %分别计算两幅图像相应分解层上对应局部区域的“能量”
E2=conv2(x2.^2,w,'same');
M=2*conv2(x1.*x2,w,'same')./(E1+E2);%计算两幅图像对应局部区域的匹配度
T=0.7;                              %定义匹配阈值
Wmin=1/2-1/2*((1-M)/(1-T));
Wmax=1-Wmin;
[m,n]=size(M);

for i=1:m
    for j=1:n
        if M(i,j)<T                %如果匹配度小于匹配阈值，说明两幅图像对应局部区域能量距离较远；
            if E1(i,j)>=E2(i,j)    %那么就直接选取区域能量较大的小波系数
                y(i,j)=x1(i,j);
            else
                y(i,j)=x2(i,j);
            end
        else                       %如果匹配度大于匹配阈值，说明两幅图像对应局部区域能量比较接近；
            if E1(i,j)>=E2(i,j)    %那么就采用加权的融合算法
                y(i,j)=Wmax(i,j)*x1(i,j)+Wmin(i,j)*x2(i,j);
            else
                y(i,j)=Wmin(i,j)*x1(i,j)+Wmax(i,j)*x2(i,j);
            end
        end
    end
end


%对比度金字塔融合函数   
function Y = fuse_con(M1, M2, zt, ap, mp)
%Y = fuse_con(M1, M2, zt, ap, mp) image fusion with contrast pyramid
%
%    M1 - input image A
%    M2 - input image B
%    zt - maximum decomposition level
%    ap - coefficient selection highpass (see selc.m) 
%    mp - coefficient selection base image (see selb.m) 
%
%    Y  - fused image   
 
%    (Oliver Rockinger 16.08.99)
 
% check inputs 
[z1,s1] = size(M1);
[z2,s2] = size(M2);
if (z1 ~= z2) || (s1 ~= s2)
  error('Input images are not of same size');
end;
 
% define filter 
w  = [1 4 6 4 1] / 16;
 
% define eps
eps = 1e-6;
 
% cells for selected images
E = cell(1,zt);
zl = zeros(1,zt);
sl = zeros(1,zt);
 
% loop over decomposition depth -> analysis
for i1 = 1:zt 
  % calculate and store actual image size 
  [z,s]  = size(M1); 
  zl(i1) = z; sl(i1)  = s; % 待计算金字塔图像的尺寸，下面判断是不是偶数，否则进行扩充1行或者1列
   
  % check if image expansion necessary 
  if (floor(z/2) ~= z/2), ew(1) = 1; else  ew(1) = 0; end;
  if (floor(s/2) ~= s/2), ew(2) = 1; else  ew(2) = 0; end;
 
  % perform expansion if necessary
  if (any(ew))
  M1 = adb(M1,ew); % adds rows and/or columns by duplication on the lower resp. right side of the input matrix
  M2 = adb(M2,ew);
  end;
   
  % perform filtering 
  G1 = conv2(conv2(es2(M1,2), w, 'valid'),w', 'valid');
  G2 = conv2(conv2(es2(M2,2), w, 'valid'),w', 'valid');
  
  % decimate, undecimate and interpolate 
  M1T = conv2(conv2(es2(undec2(dec2(G1)), 2), 2*w, 'valid'),2*w', 'valid'); % es2:symmetric extension of a matrix on all borders
  M2T = conv2(conv2(es2(undec2(dec2(G2)), 2), 2*w, 'valid'),2*w', 'valid');
  
  % select coefficients and store them
  E(i1) = {selc(M1./(M1T+eps)-1, M2./(M2T+eps)-1, ap)};
   
  % decimate 
  M1 = dec2(G1);
  M2 = dec2(G2);
end;
 
% select base coefficients of last decompostion stage
M1 = selb(M1,M2,mp);
 
% loop over decomposition depth -> synthesis
for i1 = zt:-1:1
  % undecimate and interpolate 
  M1T = conv2(conv2(es2(undec2(M1), 2), 2*w, 'valid'), 2*w', 'valid');
  % add coefficients
  M1  = (M1T+eps) .* (E{i1}+1);
  % select valid image region 
  M1 = M1(1:zl(i1),1:sl(i1));
end;
 
% copy image
Y = M1;

%lapulas function

function Y = fuse_lap(M1, M2, zt, ap, mp) 
%Y = fuse_lap(M1, M2, zt, ap, mp) image fusion with laplacian pyramid 
% 
% M1 - input image A 
% M2 - input image B 
% zt - maximum decomposition level 
% ap - coefficient selection highpass (see selc.m) 
% mp - coefficient selection base image (see selb.m) 
% 
% Y - fused image 

% (Oliver Rockinger 16.08.99) 

% check inputs 
[z1 s1] = size(M1); 
[z2 s2] = size(M2); 
if (z1 ~= z2) | (s1 ~= s2) 
error('Input images are not of same size'); 
end; 

% define filter 
w = [1 4 6 4 1] / 16; 

% cells for selected images 
E = cell(1,zt); 

% loop over decomposition depth -> analysis 
for i1 = 1:zt 
% calculate and store actual image size 
[z s] = size(M1); 
zl(i1) = z; sl(i1) = s; 

% check if image expansion necessary 
if (floor(z/2) ~= z/2), ew(1) = 1; else, ew(1) = 0; end; 
if (floor(s/2) ~= s/2), ew(2) = 1; else, ew(2) = 0; end; 

% perform expansion if necessary 
if (any(ew)) 
M1 = adb(M1,ew); 
M2 = adb(M2,ew); 
end;	 

% perform filtering 
G1 = conv2(conv2(es2(M1,2), w, 'valid'),w', 'valid'); 
G2 = conv2(conv2(es2(M2,2), w, 'valid'),w', 'valid'); 

% decimate, undecimate and interpolate 
M1T = conv2(conv2(es2(undec2(dec2(G1)), 2), 2*w, 'valid'),2*w', 'valid'); 
M2T = conv2(conv2(es2(undec2(dec2(G2)), 2), 2*w, 'valid'),2*w', 'valid'); 

% select coefficients and store them 
E(i1) = {selc(M1-M1T, M2-M2T, ap)}; 

% decimate 
M1 = dec2(G1); 
M2 = dec2(G2); 
end; 

% select base coefficients of last decompostion stage 
M1 = selb(M1,M2,mp); 

% loop over decomposition depth -> synthesis 
for i1 = zt:-1:1 
% undecimate and interpolate 
M1T = conv2(conv2(es2(undec2(M1), 2), 2*w, 'valid'), 2*w', 'valid'); 
% add coefficients 
M1 = M1T + E{i1}; 
% select valid image region 
M1 = M1(1:zl(i1),1:sl(i1)); 
end; 

% copy image 
Y = M1; 


%图像扩展函数
function Y=adb(X,bd)
[z s]=size(X);
Y=zeros(z+bd(1),s+bd(2));
Y(1:z,1:s)=X;
if(bd(1)>0)
  Y(z+1:z+bd(1),1:s)=X(z-1:-1:z-bd(1),1:s);
end;
if(bd(2)>0)
  Y(1:z,s+1:s+bd(2))=X(1:z,s-1:-1:s-bd(2));
end;
if(bd(1)>0&bd(2)>0)
  Y(z+1:z+bd(1),s+1:s+bd(2))=X(z-1:-1:z-bd(1),s-1:-1:s-bd(2));
end;

function Y = es2(X, n) 
%Y = ES2(X, n) symmetric extension of a matrix on all borders 
% 
% X - input matrix 
% n - number of rows/columns to extend 
% 
% Y - extended matrix 

% (Oliver Rockinger 16.08.99) 

[z s] = size(X); 
Y = zeros(z+2*n, s+2*n); 
Y(n+1:n+z,n:-1:1) = X(:,2:1:n+1); 
Y(n+1:n+z,n+1:1:n+s) = X; 
Y(n+1:n+z,n+s+1:1:s+2*n) = X(:,s-1:-1:s-n); 
Y(n:-1:1,n+1:s+n) = X(2:1:n+1,:); 
Y(n+z+1:1:z+2*n,n+1:s+n) = X(z-1:-1:z-n,:); 


function Y = dec2(X); 
%Y = dec2(X) downsampling of a matrix by 2 
% 
% X - input matrix 
% 
% Y - output matrix 

% (Oliver Rockinger 16.08.99) 

[a b] = size(X); 
Y = X(1:2:a, 1:2:b); 



function Y = undec2(X) 
%Y = undec2(X) upsampling of a matrix by 2 
% 
% X - input matrix 
% 
% Y - output matrix 

% (Oliver Rockinger 16.08.99) 

[z s] = size(X); 
Y = zeros(2*z, 2*s); 

Y(1:2:2*z,1:2:2*s) = X; 




function Y = selb(M1, M2, mp) 
%Y = selb(M1, M2, mp) coefficient selection for base image 
% 
% M1 - coefficients A 
% M2 - coefficients B 
% mp - switch for selection type 
% mp == 1: select A 
% mp == 2: select B 
% mp == 3: average A and B 
% 
% Y - combined coefficients 

% (Oliver Rockinger 16.08.99) 

switch (mp) 
case 1, Y = M1; 
case 2, Y = M2; 
case 3, Y = (M1 + M2) / 2; 
otherwise, error('unknown option'); 
end; 



function Y = selc(M1, M2, ap) 
%Y = selc(M1, M2, ap) coefficinet selection for highpass components 
% 
% M1 - coefficients A 
% M2 - coefficients B 
% mp - switch for selection type 
% mp == 1: choose max(abs) 
% mp == 2: salience / match measure with threshold == .75 (as proposed by Burt et al) 
% mp == 3: choose max with consistency check (as proposed by Li et al) 
% mp == 4: simple choose max 
% 
% Y - combined coefficients 

% (Oliver Rockinger 16.08.99) 

% check inputs 
[z1 s1] = size(M1); 
[z2 s2] = size(M2); 
if (z1 ~= z2) | (s1 ~= s2) 
error('Input images are not of same size'); 
end; 

% switch to method 
switch(ap(1)) 
case 1, 
% choose max(abs) 
mm = (abs(M1)) > (abs(M2)); 
Y = (mm.*M1) + ((~mm).*M2); 

case 2, 
% Burts method 
um = ap(2); th = .75; 
% compute salience 
S1 = conv2(es2(M1.*M1, floor(um/2)), ones(um), 'valid'); 
S2 = conv2(es2(M2.*M2, floor(um/2)), ones(um), 'valid'); 
% compute match 
MA = conv2(es2(M1.*M2, floor(um/2)), ones(um), 'valid'); 
MA = 2 * MA ./ (S1 + S2 + eps); 
% selection 
m1 = MA > th; m2 = S1 > S2; 
w1 = (0.5 - 0.5*(1-MA) / (1-th)); 
Y = (~m1) .* ((m2.*M1) + ((~m2).*M2)); 
Y = Y + (m1 .* ((m2.*M1.*(1-w1))+((m2).*M2.*w1) + ((~m2).*M2.*(1-w1))+((~m2).*M1.*w1))); 

case 3,	
% Lis method 
um = ap(2); 
% first step 
A1 = ordfilt2(abs(es2(M1, floor(um/2))), um*um, ones(um)); 
A2 = ordfilt2(abs(es2(M2, floor(um/2))), um*um, ones(um)); 
% second step 
mm = (conv2((A1 > A2), ones(um), 'valid')) > floor(um*um/2); 
Y = (mm.*M1) + ((~mm).*M2); 

case 4, 
% simple choose max 
mm = M1 > M2; 
Y = (mm.*M1) + ((~mm).*M2); 

otherwise, 
error('unkown option'); 
end;

% --- Executes just before lapulassp is made visible.
function lapulassp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lapulassp (see VARARGIN)

% Choose default command line output for lapulassp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lapulassp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lapulassp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)%选取视频1
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%str2='E:\SP3\';
%str='E:\SP1\';   %str中存了SP1的图像
%str1='E:\SP2\';  %str1中存了SP2的图像
[filename,pathname] =uigetfile({'*.avi';'*.mp4';'*.*'},'打开视频');
%%字符串拼接 拼装路径 以上面例子说所述 此时 srt=F:\data\1.jpg
global str4;
str4=[pathname filename];
%%打开图像
%video_file1='str4';
video1=VideoReader(str4);
frame_number1=video1.NumberOfFrames;

%% 读取视频2
%video_file='B.avi';
%video=VideoReader(video_file);
%frame_number=video.NumberOfFrames;
%% 分离图片1到SP1，分离图片2到SP2
for i=1:1
    image_name1=strcat('E:\SP1\',num2str(i));
    image_name1=strcat(image_name1,'.jpg');
    I=read(video1,i);  %读出图片
    I=imresize(I,[196,228]);
    
    imwrite(I,image_name1,'jpg');                   %写图片
    
axes(handles.axes1);
imshow(I); 
end


% --- Executes on button press in pushbutton2.%选取视频2
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%str2='E:\SP3\';
%str='E:\SP1\';   %str中存了SP1的图像
%str1='E:\SP2\';  %str1中存了SP2的图像
%video_file1='A.avi';
%video1=VideoReader(video_file1);
%frame_number1=video1.NumberOfFrames;
[filename,pathname] =uigetfile({'*.avi';'*.mp4';'*.*'},'打开视频');
%%字符串拼接 拼装路径 以上面例子说所述 此时 srt=F:\data\1.jpg
global str5;
str5=[pathname filename];
%% 读取视频2
%video_file='str5.avi';
video=VideoReader(str5);
frame_number=video.NumberOfFrames;
%% 分离图片1到SP1，分离图片2到SP2
for i=1:1
    image_name=strcat('E:\SP2\',num2str(i));
    image_name=strcat(image_name,'.jpg');
    T=read(video,i);  %读出图片
    T=imresize(T,[196,228]);                       %把图像大小确定
    
    imwrite(T,image_name,'jpg');                   %写图片
    
axes(handles.axes2);
imshow(T); 
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)%拉普拉斯金字塔
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str2='E:\SP3\';
str='E:\SP1\';   %str中存了SP1的图像
str1='E:\SP2\';  %str1中存了SP2的图像
global str4;
global str5;
%video_file1='str4';
video1=VideoReader(str4);
frame_number1=video1.NumberOfFrames;
%% 读取视频2
%video_file='str5';
video=VideoReader(str5);
frame_number=video.NumberOfFrames;
%% 分离图片1到SP1，分离图片2到SP2
for i=1:frame_number1
    image_name1=strcat('E:\SP1\',num2str(i));
    image_name1=strcat(image_name1,'.jpg');
    I=read(video1,i);  %读出图片
    I=imresize(I,[196,228]);
    
    imwrite(I,image_name1,'jpg');                   %写图片
    I=[];
    image_name=strcat('E:\SP2\',num2str(i));
    image_name=strcat(image_name,'.jpg');
    T=read(video,i);  %读出图片
    T=imresize(T,[196,228]);                       %把图像大小确定
    
    imwrite(T,image_name,'jpg');                   %写图片
    T=[];
I=imread([str,num2str(i),'.jpg']); %依次读取每一幅图像
T=imread([str1,num2str(i),'.jpg']);
im=I;
[m,n]=size(im);
%im=im2double(im);
%I=im2double(I);
%T=im2double(T);
      I = double(rgb2gray(I))/255; 
      T = double(rgb2gray(T))/255; 
%普拉斯金塔变换参数 
mp = 1;zt =4; cf =1;ar = 1; cc = [cf ar]; 

Y_lap = fuse_lap(I,T,zt,cc,mp); 
imwrite(im,[str2,num2str(i),'.jpg']);
axes(handles.axes1);
imshow(I);
axes(handles.axes2);
imshow(T);
axes(handles.axes3);
imshow(Y_lap);
end
clear all;


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)%加权平均
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str2='E:\SP3\';
str='E:\SP1\';   %str中存了SP1的图像
str1='E:\SP2\';  %str1中存了SP2的图像
global str4;
video1=VideoReader(str4);
frame_number1=video1.NumberOfFrames;
%% 读取视频2
global str5;
video=VideoReader(str5);
frame_number=video.NumberOfFrames;
%% 分离图片1到SP1，分离图片2到SP2
for i=1:frame_number1
    image_name1=strcat('E:\SP1\',num2str(i));
    image_name1=strcat(image_name1,'.jpg');
    I=read(video1,i);  %读出图片
    I=imresize(I,[196,228]);
    
    imwrite(I,image_name1,'jpg');                   %写图片
    I=[];
    image_name=strcat('E:\SP2\',num2str(i));
    image_name=strcat(image_name,'.jpg');
    T=read(video,i);  %读出图片
    T=imresize(T,[196,228]);                       %把图像大小确定
    
    imwrite(T,image_name,'jpg');                   %写图片
    T=[];
I=imread([str,num2str(i),'.jpg']); %依次读取每一幅图像
T=imread([str1,num2str(i),'.jpg']);
im=I;
[m,n]=size(im);
im=im2double(im);
I=im2double(I);
T=im2double(T);
for p=1:m
    for j=1:n
        im(p,j)=I(p,j)+T(p,j);
    end
end
imwrite(im,[str2,num2str(i),'.jpg']);
axes(handles.axes1);
imshow(I);
axes(handles.axes2);
imshow(T);
axes(handles.axes3);
imshow(im); 
end
clear all;


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)%  wt小波变换
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str2='E:\SP3\';
str='E:\SP1\';   %str中存了SP1的图像
str1='E:\SP2\';  %str1中存了SP2的图像
global str4;
global str5;
%video_file1='str4';
video1=VideoReader(str4);
frame_number1=video1.NumberOfFrames;
%% 读取视频2
%video_file='str5';
video=VideoReader(str5);
frame_number=video.NumberOfFrames;
%% 分离图片1到SP1，分离图片2到SP2
for i=1:frame_number1
    image_name1=strcat('E:\SP1\',num2str(i));
    image_name1=strcat(image_name1,'.jpg');
    I=read(video1,i);  %读出图片
    I=imresize(I,[196,228]);
    
    imwrite(I,image_name1,'jpg');                   %写图片
    I=[];
    image_name=strcat('E:\SP2\',num2str(i));
    image_name=strcat(image_name,'.jpg');
    T=read(video,i);  %读出图片
    T=imresize(T,[196,228]);                       %把图像大小确定
    
    imwrite(T,image_name,'jpg');                   %写图片
    T=[];
I=imread([str,num2str(i),'.jpg']); %依次读取每一幅图像
T=imread([str1,num2str(i),'.jpg']);
im=I;
[m,n]=size(im);
%im=im2double(im);
%I=im2double(I);
%T=im2double(T);
      I = double(rgb2gray(I))/255;         %小波需要这个吗
      T = double(rgb2gray(T))/255;         %小波需要这个吗

Y_wt=wtfusion(I,T,2,'db1'); 
imwrite(im,[str2,num2str(i),'.jpg']);
axes(handles.axes1);
imshow(I);
axes(handles.axes2);
imshow(T);
axes(handles.axes3);
imshow(Y_wt);
end
clear all;


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)%RP金字塔
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str2='E:\SP3\';
str='E:\SP1\';   %str中存了SP1的图像
str1='E:\SP2\';  %str1中存了SP2的图像
global str4;
global str5;
%video_file1='str4';
video1=VideoReader(str4);
frame_number1=video1.NumberOfFrames;
%% 读取视频2
%video_file='str5';
video=VideoReader(str5);
frame_number=video.NumberOfFrames;
%% 分离图片1到SP1，分离图片2到SP2
for i=1:frame_number1
    image_name1=strcat('E:\SP1\',num2str(i));
    image_name1=strcat(image_name1,'.jpg');
    I=read(video1,i);  %读出图片
    I=imresize(I,[196,228]);
    
    imwrite(I,image_name1,'jpg');                   %写图片
    I=[];
    image_name=strcat('E:\SP2\',num2str(i));
    image_name=strcat(image_name,'.jpg');
    T=read(video,i);  %读出图片
    T=imresize(T,[196,228]);                       %把图像大小确定
    
    imwrite(T,image_name,'jpg');                   %写图片
    T=[];
I=imread([str,num2str(i),'.jpg']); %依次读取每一幅图像
T=imread([str1,num2str(i),'.jpg']);
im=I;
[m,n]=size(im);
%im=im2double(im);
%I=im2double(I);
%T=im2double(T);
      I = double(rgb2gray(I))/255;      %   RP金字塔需要这个吗
      T = double(rgb2gray(T))/255;      %   RP金字塔需要这个吗
%	RP金塔变换参数 
mp = 1;zt =4; cf =1;ar = 1; cc = [cf ar]; 

Y_RP= RP_fuse(I,T,zt,cc,mp); 
imwrite(im,[str2,num2str(i),'.jpg']);
axes(handles.axes1);
imshow(I);
axes(handles.axes2);
imshow(T);
axes(handles.axes3);
imshow(Y_RP);
end
clear all;

% --- Executes on button press in pushbutton3.
function pushbutton7_Callback(hObject, eventdata, handles)%对比金字塔
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str2='E:\SP3\';
str='E:\SP1\';   %str中存了SP1的图像
str1='E:\SP2\';  %str1中存了SP2的图像
global str4;
global str5;
%video_file1='str4';
video1=VideoReader(str4);
frame_number1=video1.NumberOfFrames;
%% 读取视频2
%video_file='str5';
video=VideoReader(str5);
frame_number=video.NumberOfFrames;
%% 分离图片1到SP1，分离图片2到SP2
for i=1:frame_number1
    image_name1=strcat('E:\SP1\',num2str(i));
    image_name1=strcat(image_name1,'.jpg');
    I=read(video1,i);  %读出图片
    I=imresize(I,[196,228]);
    
    imwrite(I,image_name1,'jpg');                   %写图片
    I=[];
    image_name=strcat('E:\SP2\',num2str(i));
    image_name=strcat(image_name,'.jpg');
    T=read(video,i);  %读出图片
    T=imresize(T,[196,228]);                       %把图像大小确定
    
    imwrite(T,image_name,'jpg');                   %写图片
    T=[];
I=imread([str,num2str(i),'.jpg']); %依次读取每一幅图像
T=imread([str1,num2str(i),'.jpg']);
im=I;
[m,n]=size(im);
%im=im2double(im);
%I=im2double(I);
%T=im2double(T);
      I = double(rgb2gray(I))/255;     %为什么是这个
      T = double(rgb2gray(T))/255;     %为什么是这个
%对比金塔变换参数 
mp = 1;zt =4; cf =1;ar = 1; cc = [cf ar]; 

Y_con=fuse_con(I,T,zt,cc,mp); 
imwrite(im,[str2,num2str(i),'.jpg']);
axes(handles.axes1);
imshow(I);
axes(handles.axes2);
imshow(T);
axes(handles.axes3);
imshow(Y_con);
end
clear all;
