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


%RP������

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

%С���任wt

function y= wtfusion(x1,x2,N,wname)

%�������ܣ�
%     ����x=wtfusion(x1,x2,N,wname)������ԭͼ��x1,x2���л���С���任��ͼ���ںϣ��õ��ںϺ��ͼ��y
%     ���Ʒ������ü�Ȩƽ�����ںϹ��򣬸�ϸ�ڷ������û�����������������ںϹ���
%���������
%     x1----����ԭͼ��1
%     x2----����ԭͼ��2
%     N----С���ֽ�Ĳ���
%     wname----С��������
%���������
%     y----ԭͼ���ںϺ�õ���ͼ��
%-----------------------------------------------------------------%

x1=double(x1);                   %��uint8��ͼ����������ת����double�ͽ������ݴ���
x2=double(x2);

 %��ԭͼ��x1,x2�ֱ����N��С���ֽ⣬wnameΪС����������
 %CΪ����ֽ�ϵ��,SΪ����ֽ�ϵ������,Ҳ���Ǵ�С. 
 %C�Ľṹ:c=[A(N)|H(N)|V(N)|D(N)|H(N-1)|V(N-1)|D(N-1)|H(N-2)|V(N-2)|D(N-2)|...|H(1)|V(1)|D(1)]
 %A(N)�����N���Ƶϵ��,H(N)|V(N)|D(N)�����N���Ƶϵ��,�ֱ���ˮƽ,��ֱ,�ԽǸ�Ƶ
 %
 %S��һ��(N+2)�С�2�еĽṹ�Ǵ������ֽ�ϵ�����ȵ�,����һ����A(N)�ĳ��ȣ���ʵ��A(N)��ԭ�����������������,
 %�ڶ�����H(N)|V(N)|D(N)|�ĳ���,��������H(N-1)|V(N-1)|D(N-1)�ĳ���,
 %�����ڶ�����H(1)|V(1)|D(1)����,���һ����X�ĳ���(��С)

[C1,S1]=wavedec2(x1,N,wname);    % �������ܣ�ʵ��ͼ���С���ֽ⡣wname��С���������ơ�N�Ƿֽ�Ĳ���.
[C2,S2]=wavedec2(x2,N,wname);    

A1=appcoef2(C1,S1,wname,N);            %�������ܣ���ɢС���任��Ƶ����ϵ����ȡ
A2=appcoef2(C2,S2,wname,N);             %��ȡ��С���ֽ�Ľ��Ʒ�������Ƶֻ��һ��������A(N)������Ķ��Ǹ�Ƶ��
A=0.5*A1+0.5*A2;                       %��Ƶ�������ںϹ�����ü�Ȩƽ���ķ���

%figure
%imshow(uint8(A))

%����matlab�н��Ʒ�����ϸ�ڷ����Ĵ洢��ʽ�����ںϺ�Ľ��Ʒ�����ϸ�ڷ���ת����������Ȼ���������C��
%��������Ϊ�˷����ع�ԭͼ��

a=reshape(A,1,S1(1,1)*S1(1,2));        %��Aת����������.
C=a;

for i=N:-1:1                           %ѭ���ӵ�N�㵽��1��    
    [H1,V1,D1]=detcoef2('all',C1,S1,i);       %С���任��Ƶ������ȡ
    [H2,V2,D2]=detcoef2('all',C2,S2,i);        %��ȡ��С���ֽ�ĸ���ϸ�ڷ���
    H=f(H1,H2);                           
    V=f(V1,V2);    %��Ƶ�������ںϹ��򣺲��û�����������������ںϹ��򡣿���ѡ��ͬ���ںϹ���
    D=f(D1,D2);
    h=reshape(H,1,S1(N+2-i,1)*S1(N+2-i,2));%�ֱ��ںϺ��ϸ�ڷ���ת����������������������C��
    v=reshape(V,1,S1(N+2-i,1)*S1(N+2-i,2));
    d=reshape(D,1,S1(N+2-i,1)*S1(N+2-i,2));
    C=[C,h,v,d];
end

S=S1;
y=waverec2(C,S,wname);      %�ع�ԭͼ��
%figure(1);imshow(uint8(y));title('����С���任���ں�ͼ��')


%С���任f�Ӻ���
function y = f(x1 , x2)

%�������ܣ�
%       y=f(x1,x2)������ԭͼ��x1��x2������������������ںϹ�������ں�,�õ��ںϺ��ͼ��y
%       ���ȼ�������ͼ���ƥ��ȣ���ƥ��ȴ�����ֵ��˵������ͼ���Ӧ�ֲ������Ͻӽ���
%       ��˲��ü�Ȩƽ�����ںϷ�������ƥ���С����ֵ��˵������ͼ���Ӧ�ֲ��������ϴ�
%       ���ѡȡ�ֲ����������ϴ��С��ϵ����Ϊ�ں�ͼ���С��ϵ��
%���������
%      x1----����ԭͼ��1
%      x2----����ԭͼ��2
%���������
%      y----�ںϺ��ͼ��
%------------------------------------------------------------%

w=1/16*[1 2 1;2 4 2;1 2 1];   %Ȩϵ��
E1=conv2(x1.^2,w,'same');     %�ֱ��������ͼ����Ӧ�ֽ���϶�Ӧ�ֲ�����ġ�������
E2=conv2(x2.^2,w,'same');
M=2*conv2(x1.*x2,w,'same')./(E1+E2);%��������ͼ���Ӧ�ֲ������ƥ���
T=0.7;                              %����ƥ����ֵ
Wmin=1/2-1/2*((1-M)/(1-T));
Wmax=1-Wmin;
[m,n]=size(M);

for i=1:m
    for j=1:n
        if M(i,j)<T                %���ƥ���С��ƥ����ֵ��˵������ͼ���Ӧ�ֲ��������������Զ��
            if E1(i,j)>=E2(i,j)    %��ô��ֱ��ѡȡ���������ϴ��С��ϵ��
                y(i,j)=x1(i,j);
            else
                y(i,j)=x2(i,j);
            end
        else                       %���ƥ��ȴ���ƥ����ֵ��˵������ͼ���Ӧ�ֲ����������ȽϽӽ���
            if E1(i,j)>=E2(i,j)    %��ô�Ͳ��ü�Ȩ���ں��㷨
                y(i,j)=Wmax(i,j)*x1(i,j)+Wmin(i,j)*x2(i,j);
            else
                y(i,j)=Wmin(i,j)*x1(i,j)+Wmax(i,j)*x2(i,j);
            end
        end
    end
end


%�ԱȶȽ������ںϺ���   
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
  zl(i1) = z; sl(i1)  = s; % �����������ͼ��ĳߴ磬�����ж��ǲ���ż���������������1�л���1��
   
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


%ͼ����չ����
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
function pushbutton1_Callback(hObject, eventdata, handles)%ѡȡ��Ƶ1
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%str2='E:\SP3\';
%str='E:\SP1\';   %str�д���SP1��ͼ��
%str1='E:\SP2\';  %str1�д���SP2��ͼ��
[filename,pathname] =uigetfile({'*.avi';'*.mp4';'*.*'},'����Ƶ');
%%�ַ���ƴ�� ƴװ·�� ����������˵���� ��ʱ srt=F:\data\1.jpg
global str4;
str4=[pathname filename];
%%��ͼ��
%video_file1='str4';
video1=VideoReader(str4);
frame_number1=video1.NumberOfFrames;

%% ��ȡ��Ƶ2
%video_file='B.avi';
%video=VideoReader(video_file);
%frame_number=video.NumberOfFrames;
%% ����ͼƬ1��SP1������ͼƬ2��SP2
for i=1:1
    image_name1=strcat('E:\SP1\',num2str(i));
    image_name1=strcat(image_name1,'.jpg');
    I=read(video1,i);  %����ͼƬ
    I=imresize(I,[196,228]);
    
    imwrite(I,image_name1,'jpg');                   %дͼƬ
    
axes(handles.axes1);
imshow(I); 
end


% --- Executes on button press in pushbutton2.%ѡȡ��Ƶ2
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%str2='E:\SP3\';
%str='E:\SP1\';   %str�д���SP1��ͼ��
%str1='E:\SP2\';  %str1�д���SP2��ͼ��
%video_file1='A.avi';
%video1=VideoReader(video_file1);
%frame_number1=video1.NumberOfFrames;
[filename,pathname] =uigetfile({'*.avi';'*.mp4';'*.*'},'����Ƶ');
%%�ַ���ƴ�� ƴװ·�� ����������˵���� ��ʱ srt=F:\data\1.jpg
global str5;
str5=[pathname filename];
%% ��ȡ��Ƶ2
%video_file='str5.avi';
video=VideoReader(str5);
frame_number=video.NumberOfFrames;
%% ����ͼƬ1��SP1������ͼƬ2��SP2
for i=1:1
    image_name=strcat('E:\SP2\',num2str(i));
    image_name=strcat(image_name,'.jpg');
    T=read(video,i);  %����ͼƬ
    T=imresize(T,[196,228]);                       %��ͼ���Сȷ��
    
    imwrite(T,image_name,'jpg');                   %дͼƬ
    
axes(handles.axes2);
imshow(T); 
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)%������˹������
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str2='E:\SP3\';
str='E:\SP1\';   %str�д���SP1��ͼ��
str1='E:\SP2\';  %str1�д���SP2��ͼ��
global str4;
global str5;
%video_file1='str4';
video1=VideoReader(str4);
frame_number1=video1.NumberOfFrames;
%% ��ȡ��Ƶ2
%video_file='str5';
video=VideoReader(str5);
frame_number=video.NumberOfFrames;
%% ����ͼƬ1��SP1������ͼƬ2��SP2
for i=1:frame_number1
    image_name1=strcat('E:\SP1\',num2str(i));
    image_name1=strcat(image_name1,'.jpg');
    I=read(video1,i);  %����ͼƬ
    I=imresize(I,[196,228]);
    
    imwrite(I,image_name1,'jpg');                   %дͼƬ
    I=[];
    image_name=strcat('E:\SP2\',num2str(i));
    image_name=strcat(image_name,'.jpg');
    T=read(video,i);  %����ͼƬ
    T=imresize(T,[196,228]);                       %��ͼ���Сȷ��
    
    imwrite(T,image_name,'jpg');                   %дͼƬ
    T=[];
I=imread([str,num2str(i),'.jpg']); %���ζ�ȡÿһ��ͼ��
T=imread([str1,num2str(i),'.jpg']);
im=I;
[m,n]=size(im);
%im=im2double(im);
%I=im2double(I);
%T=im2double(T);
      I = double(rgb2gray(I))/255; 
      T = double(rgb2gray(T))/255; 
%����˹�����任���� 
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
function pushbutton4_Callback(hObject, eventdata, handles)%��Ȩƽ��
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str2='E:\SP3\';
str='E:\SP1\';   %str�д���SP1��ͼ��
str1='E:\SP2\';  %str1�д���SP2��ͼ��
global str4;
video1=VideoReader(str4);
frame_number1=video1.NumberOfFrames;
%% ��ȡ��Ƶ2
global str5;
video=VideoReader(str5);
frame_number=video.NumberOfFrames;
%% ����ͼƬ1��SP1������ͼƬ2��SP2
for i=1:frame_number1
    image_name1=strcat('E:\SP1\',num2str(i));
    image_name1=strcat(image_name1,'.jpg');
    I=read(video1,i);  %����ͼƬ
    I=imresize(I,[196,228]);
    
    imwrite(I,image_name1,'jpg');                   %дͼƬ
    I=[];
    image_name=strcat('E:\SP2\',num2str(i));
    image_name=strcat(image_name,'.jpg');
    T=read(video,i);  %����ͼƬ
    T=imresize(T,[196,228]);                       %��ͼ���Сȷ��
    
    imwrite(T,image_name,'jpg');                   %дͼƬ
    T=[];
I=imread([str,num2str(i),'.jpg']); %���ζ�ȡÿһ��ͼ��
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
function pushbutton5_Callback(hObject, eventdata, handles)%  wtС���任
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str2='E:\SP3\';
str='E:\SP1\';   %str�д���SP1��ͼ��
str1='E:\SP2\';  %str1�д���SP2��ͼ��
global str4;
global str5;
%video_file1='str4';
video1=VideoReader(str4);
frame_number1=video1.NumberOfFrames;
%% ��ȡ��Ƶ2
%video_file='str5';
video=VideoReader(str5);
frame_number=video.NumberOfFrames;
%% ����ͼƬ1��SP1������ͼƬ2��SP2
for i=1:frame_number1
    image_name1=strcat('E:\SP1\',num2str(i));
    image_name1=strcat(image_name1,'.jpg');
    I=read(video1,i);  %����ͼƬ
    I=imresize(I,[196,228]);
    
    imwrite(I,image_name1,'jpg');                   %дͼƬ
    I=[];
    image_name=strcat('E:\SP2\',num2str(i));
    image_name=strcat(image_name,'.jpg');
    T=read(video,i);  %����ͼƬ
    T=imresize(T,[196,228]);                       %��ͼ���Сȷ��
    
    imwrite(T,image_name,'jpg');                   %дͼƬ
    T=[];
I=imread([str,num2str(i),'.jpg']); %���ζ�ȡÿһ��ͼ��
T=imread([str1,num2str(i),'.jpg']);
im=I;
[m,n]=size(im);
%im=im2double(im);
%I=im2double(I);
%T=im2double(T);
      I = double(rgb2gray(I))/255;         %С����Ҫ�����
      T = double(rgb2gray(T))/255;         %С����Ҫ�����

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
function pushbutton6_Callback(hObject, eventdata, handles)%RP������
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str2='E:\SP3\';
str='E:\SP1\';   %str�д���SP1��ͼ��
str1='E:\SP2\';  %str1�д���SP2��ͼ��
global str4;
global str5;
%video_file1='str4';
video1=VideoReader(str4);
frame_number1=video1.NumberOfFrames;
%% ��ȡ��Ƶ2
%video_file='str5';
video=VideoReader(str5);
frame_number=video.NumberOfFrames;
%% ����ͼƬ1��SP1������ͼƬ2��SP2
for i=1:frame_number1
    image_name1=strcat('E:\SP1\',num2str(i));
    image_name1=strcat(image_name1,'.jpg');
    I=read(video1,i);  %����ͼƬ
    I=imresize(I,[196,228]);
    
    imwrite(I,image_name1,'jpg');                   %дͼƬ
    I=[];
    image_name=strcat('E:\SP2\',num2str(i));
    image_name=strcat(image_name,'.jpg');
    T=read(video,i);  %����ͼƬ
    T=imresize(T,[196,228]);                       %��ͼ���Сȷ��
    
    imwrite(T,image_name,'jpg');                   %дͼƬ
    T=[];
I=imread([str,num2str(i),'.jpg']); %���ζ�ȡÿһ��ͼ��
T=imread([str1,num2str(i),'.jpg']);
im=I;
[m,n]=size(im);
%im=im2double(im);
%I=im2double(I);
%T=im2double(T);
      I = double(rgb2gray(I))/255;      %   RP��������Ҫ�����
      T = double(rgb2gray(T))/255;      %   RP��������Ҫ�����
%	RP�����任���� 
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
function pushbutton7_Callback(hObject, eventdata, handles)%�ԱȽ�����
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str2='E:\SP3\';
str='E:\SP1\';   %str�д���SP1��ͼ��
str1='E:\SP2\';  %str1�д���SP2��ͼ��
global str4;
global str5;
%video_file1='str4';
video1=VideoReader(str4);
frame_number1=video1.NumberOfFrames;
%% ��ȡ��Ƶ2
%video_file='str5';
video=VideoReader(str5);
frame_number=video.NumberOfFrames;
%% ����ͼƬ1��SP1������ͼƬ2��SP2
for i=1:frame_number1
    image_name1=strcat('E:\SP1\',num2str(i));
    image_name1=strcat(image_name1,'.jpg');
    I=read(video1,i);  %����ͼƬ
    I=imresize(I,[196,228]);
    
    imwrite(I,image_name1,'jpg');                   %дͼƬ
    I=[];
    image_name=strcat('E:\SP2\',num2str(i));
    image_name=strcat(image_name,'.jpg');
    T=read(video,i);  %����ͼƬ
    T=imresize(T,[196,228]);                       %��ͼ���Сȷ��
    
    imwrite(T,image_name,'jpg');                   %дͼƬ
    T=[];
I=imread([str,num2str(i),'.jpg']); %���ζ�ȡÿһ��ͼ��
T=imread([str1,num2str(i),'.jpg']);
im=I;
[m,n]=size(im);
%im=im2double(im);
%I=im2double(I);
%T=im2double(T);
      I = double(rgb2gray(I))/255;     %Ϊʲô�����
      T = double(rgb2gray(T))/255;     %Ϊʲô�����
%�ԱȽ����任���� 
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
