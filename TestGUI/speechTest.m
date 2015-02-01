function varargout = speechTest(varargin)
% SPEECHTEST MATLAB code for speechTest.fig
%      SPEECHTEST, by itself, creates a new SPEECHTEST or raises the existing
%      singleton*.
%
%      H = SPEECHTEST returns the handle to a new SPEECHTEST or the handle to
%      the existing singleton*.
%
%      SPEECHTEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPEECHTEST.M with the given input arguments.
%
%      SPEECHTEST('Property','Value',...) creates a new SPEECHTEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before speechTest_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to speechTest_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help speechTest

% Last Modified by GUIDE v2.5 17-Mar-2014 15:39:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @speechTest_OpeningFcn, ...
                   'gui_OutputFcn',  @speechTest_OutputFcn, ...
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


% --- Executes just before speechTest is made visible.
function speechTest_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to speechTest (see VARARGIN)

% Choose default command line output for speechTest
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes speechTest wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = speechTest_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonStart.
function pushbuttonStart_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This function sets things up and starts the test. It will first index the
% directory to which you have saved all the different sound files processed
% through various conditions like noisy, idbm, audio, audio-video and
% video. Then it will pick out only ones that have not been degraded by
% initial processing and have no context. Then it will present the
% sentences in random order and such that no sentence is repeated twice
% during presentation. All results will be displayed on the uitable. 
% To quickly visualize the results use the evaluate button. 
% To write data to file press the write button.
% Please configure the path before pressing the start button or it won't work.
% Author: Arun.P.U.
    %testTypeval = get(handles.popupmenu2,'Value');
    testTypeval = 2; % bydefault this is type two. Type 1 was introduced 
    % with the idea in mind that we may want to check at which SNR the subject
    % falls below a certain comprehension threhsold and presenting the
    % actual test at that SNR. Howver this idea was later dropped and the SNr was fixed at -8.
    SNRval = get(handles.popupmenuSNR,'Value');
    SNR=-1*(2+SNRval*-2);
    soundPath = 'C:\Work\AV_project\code\BinaryMask\idbm\TestGUI\Old\test_CN_-8_40_59\SNR_8dB';% this is where the sound files reside
    %dirPath = [soundPath,'\SNR_',num2str(SNR),'dB'];
    dirPath = soundPath;
    sentPath = 'C:\Work\AV_project\data\sentences'; % put the sound files in a subdirectory and cponfigure this path
    sentDir = dir([sentPath,'\*.txt']);
    restdata = separateFiles(sentDir);% get those indices that have no context
    % and have not been corrupted by SSBoll79.
    %restdata = setdiff(restdata,nocontIndex);
    %restIndex=1:201;
    set(handles.pushbuttonStart,'UserData',restdata);
    
    switch testTypeval        
        case 1
            n=5;
            % Index the noisy files and save rest
            soundFilesN =dir([dirPath,'\Noisy*.wav']);
            restIndex=get(handles.pushbuttonStart,'UserData'); 
            fileIndexN=randsample(restIndex,n);
            restIndex = setdiff(restIndex,fileIndexN);
            set(handles.pushbuttonStart,'UserData',restIndex);
            % Index Idbm and save rest
            soundFilesI =dir([dirPath,'\IDBM*.wav']);
            fileIndexI=randsample(restIndex,n);
            restIndex = setdiff(restIndex,fileIndexI);
            set(handles.pushbuttonStart,'UserData',restIndex);
            n1=randperm(n);
            for j = 1:n
                set(handles.text2,'String',['Sentence ',num2str(j),' of ',num2str(n)]);
                i = n1(j);
                % play a noisy file
                if mod(i,2)==1                    
                    soundName=getfield( soundFilesN,{fileIndexN(i)},'name');
                    soundfilepath=[dirPath,'\',soundName];
                    [wavefile,fs]=wavread(soundfilepath);
                    sound(wavefile,fs);
                    %display the sentence
                    sentName=getfield( sentDir,{fileIndexN(i)},'name');
                    sentID1=strsplit(sentName,'-');
                    sentID=str2num(sentID1{1});
                    set(handles.pushbuttonEval,'UserData',sentID);
                    sentfilepath=[sentPath,'\',sentName];
                    fid=fopen(sentfilepath);
                    C=textscan(fid,'%s');
                    C1=C{1}';
                    set(handles.uitable2,'Data',C1)
                    set(handles.pushbuttonStart,'UserData','Noisy');
                    fclose(fid);
                    uiwait;
                else                
                % play an IDBM
                    soundName=getfield( soundFilesI,{fileIndexI(i)},'name');
                     soundfilepath=[dirPath,'\',soundName];
                    [wavefile,fs]=wavread(soundfilepath);
                    sound(wavefile,fs);
                    %display the sentence
                    sentName=getfield( sentDir,{fileIndexI(i)},'name');
                    sentID1=strsplit(sentName,'-');
                    sentID=str2num(sentID1{1});
                    set(handles.pushbuttonEval,'UserData',sentID);
                    sentfilepath=[sentPath,'\',sentName];
                    fid=fopen(sentfilepath);
                    C=textscan(fid,'%s');
                    C1=C{1}';
                    set(handles.uitable2,'Data',C1)
                    set(handles.pushbuttonStart,'UserData','IDBM');
                    fclose(fid);
                    uiwait;
                end               
            end            
        case 2
             % Here we index files for a certain condition and save the rest
             % so that no sentence is repeated for a second time.
            n=10;% number of senteces for each condition
            total = n*5;% total number of sentences presented
            % Index the noisy files and save rest
            soundFilesN =dir([dirPath,'\Noisy*.wav']);
            restIndex=get(handles.pushbuttonStart,'UserData'); 
            fileIndexN=randsample(restIndex,n);
            restIndex = setdiff(restIndex,fileIndexN);
            set(handles.pushbuttonStart,'UserData',restIndex);
            % Index Idbm and save rest
            soundFilesI =dir([dirPath,'\IDBM*.wav']);
            fileIndexI=randsample(restIndex,n);
            restIndex = setdiff(restIndex,fileIndexI);
            set(handles.pushbuttonStart,'UserData',restIndex);
            % Index Audio and save rest
            soundFilesA =dir([dirPath,'\Audio_*.wav']);
            fileIndexA=randsample(restIndex,n);
            restIndex = setdiff(restIndex,fileIndexA);
            set(handles.pushbuttonStart,'UserData',restIndex);
            % Index Audio Video and save rest
            soundFilesAV =dir([dirPath,'\AudioVideo*.wav']);
            fileIndexAV=randsample(restIndex,n);
            restIndex = setdiff(restIndex,fileIndexAV);
            set(handles.pushbuttonStart,'UserData',restIndex);
            % Index Audio Video and save rest
            soundFilesV =dir([dirPath,'\Video*.wav']);
            fileIndexV=randsample(restIndex,n);
            restIndex = setdiff(restIndex,fileIndexV);
            set(handles.pushbuttonStart,'UserData',restIndex); 
            n1=randperm(total);
            countn=1;% initialize counts, these keep track of 
            counti=1;% how many sound files have been played
            counta=1;% for any given condition
            countav=1;% 
            countv=1;% 
            for j = 1:total
                set(handles.text2,'String',['Sentence ',num2str(j),' of ',num2str(total)]);
                i = n1(j);
                % play a noisy file
                switch mod(i,5)
                case 1                   
                    soundName=getfield( soundFilesN,{fileIndexN(countn)},'name');
                    soundfilepath=[dirPath,'\',soundName];
                    [wavefile,fs]=wavread(soundfilepath);
                    sound(wavefile,fs);
                    %display the sentence
                    sentName=getfield( sentDir,{fileIndexN(countn)},'name');
                    sentID1=strsplit(sentName,'-');
                    sentID=str2num(sentID1{1});
                    set(handles.pushbuttonEval,'UserData',sentID);
                    sentfilepath=[sentPath,'\',sentName];
                    fid=fopen(sentfilepath);
                    C=textscan(fid,'%s');
                    C1=C{1}';
                    set(handles.uitable2,'Data',C1)
                    set(handles.pushbuttonStart,'UserData','Noisy');
                    fclose(fid);
                    uiwait;
                    countn=countn+1;
                case 2                
                % play an IDBM
                    soundName=getfield( soundFilesI,{fileIndexI(counti)},'name');
                    soundfilepath=[dirPath,'\',soundName];
                    [wavefile,fs]=wavread(soundfilepath);
                    sound(wavefile,fs);
                    %display the sentence
                    sentName=getfield( sentDir,{fileIndexI(counti)},'name');
                    sentID1=strsplit(sentName,'-');
                    sentID=str2num(sentID1{1});
                    set(handles.pushbuttonEval,'UserData',sentID);
                    sentfilepath=[sentPath,'\',sentName];
                    fid=fopen(sentfilepath);
                    C=textscan(fid,'%s');
                    C1=C{1}';
                    set(handles.uitable2,'Data',C1)
                    set(handles.pushbuttonStart,'UserData','IDBM');
                    fclose(fid);
                    uiwait;
                    counti=counti+1;
                case 3                
                % play an Audio file
                    soundName=getfield( soundFilesA,{fileIndexA(counta)},'name');
                    soundfilepath=[dirPath,'\',soundName];
                    [wavefile,fs]=wavread(soundfilepath);
                    sound(wavefile,fs);
                    %display the sentence
                    sentName=getfield( sentDir,{fileIndexA(counta)},'name');
                    sentID1=strsplit(sentName,'-');
                    sentID=str2num(sentID1{1});
                    set(handles.pushbuttonEval,'UserData',sentID);
                    sentfilepath=[sentPath,'\',sentName];
                    fid=fopen(sentfilepath);
                    C=textscan(fid,'%s');
                    C1=C{1}';
                    set(handles.uitable2,'Data',C1)
                    set(handles.pushbuttonStart,'UserData','Audio');
                    fclose(fid);
                    uiwait;
                    counta=counta+1;
                case 4  
                % play an AudioVideo file
                    soundName=getfield( soundFilesAV,{fileIndexAV(countav)},'name');
                    soundfilepath=[dirPath,'\',soundName];
                    [wavefile,fs]=wavread(soundfilepath);
                    sound(wavefile,fs);
                    %display the sentence
                    sentName=getfield( sentDir,{fileIndexAV(countav)},'name');
                    sentID1=strsplit(sentName,'-');
                    sentID=str2num(sentID1{1});
                    set(handles.pushbuttonEval,'UserData',sentID);
                    sentfilepath=[sentPath,'\',sentName];
                    fid=fopen(sentfilepath);
                    C=textscan(fid,'%s');
                    C1=C{1}';
                    set(handles.uitable2,'Data',C1)
                    set(handles.pushbuttonStart,'UserData','AudioVideo');
                    fclose(fid);
                    uiwait;
                    countav=countav+1;
                case 0 
                % play a Video file
                    soundName=getfield( soundFilesV,{fileIndexV(countv)},'name');
                    soundfilepath=[dirPath,'\',soundName];
                    [wavefile,fs]=wavread(soundfilepath);
                    sound(wavefile,fs);
                    %display the sentence
                    sentName=getfield( sentDir,{fileIndexV(countv)},'name');
                    sentID1=strsplit(sentName,'-');
                    sentID=str2num(sentID1{1});
                    set(handles.pushbuttonEval,'UserData',sentID);
                    sentfilepath=[sentPath,'\',sentName];
                    fid=fopen(sentfilepath);
                    C=textscan(fid,'%s');
                    C1=C{1}';
                    set(handles.uitable2,'Data',C1)
                    set(handles.pushbuttonStart,'UserData','Video');
                    fclose(fid);
                    uiwait;
                    countv=countv+1;
                end
                           
            end
    end


% --- Executes on button press in pushbuttonRepeat.
function pushbuttonRepeat_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRepeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonStop.
function pushbuttonStop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume; 
close all;


% --- Executes on selection change in popupmenuSNR.
function popupmenuSNR_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuSNR contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuSNR


% --- Executes during object creation, after setting all properties.
function popupmenuSNR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected cell(s) is changed in uitable2.
function uitable2_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitable2,'UserData',eventdata.Indices);
%eventdata.Indices


% --- Executes on button press in pushbuttonEval.
function pushbuttonEval_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This function evaluates the total number of words identified correctly for any sentence presented.
% And stores the result in uitable3.UserData and uitable3.Data.
% Author: Arun.P.U.
    uiresume
    data2=get(handles.uitable3,'UserData');
    l1=get(handles.uitable2,'Data');
    l2=get(handles.uitable2,'UserData');
    data1{1,1}=get(handles.pushbuttonEval,'UserData');
    data1{1,2}=size(l1,2);
    data1{1,3}=size(l2,1);
    data1{1,4}=(size(l2,1)/length(l1))*100;

    %set(handles.uitable2,'UserData',eventdata.Indices);
    % set(handles.uitable2,'UserData',C1);
    string1=[];
    for i =1:length(l1)
       string1=[string1,';',l1{i}];
    end
    string2=[];
    for i =1:size(l2,1)
       if isempty(l2)
           string2=[];
        elseif size(l2,2)==1
           string2=[string2,';',l1{l2(2)}];       
       else       
           string2=[string2,';',l1{l2(i,2)}];
       end
    end
    data1{1,5}=string1;
    data1{1,6}=string2;    
    dataType=get(handles.pushbuttonStart,'UserData');
    data1{1,7} = dataType;    
    data2=[data2;data1];
    set(handles.uitable3,'UserData',data2);
    set(handles.uitable3,'Data',data2);
   
        
    
    
% --- Executes on button press in pushbuttonWrite.
function pushbuttonWrite_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonWrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% It is useful to store the data collected in the test and displayed and
% stored inuitable to be directly written to  a csv file. I prefer
% importing this csv file to R and plotting the results.
% Author: Arun.P.U.
writeFile=1;
    if writeFile        
        resultsPath='C:\Work\AV_project\code\BinaryMask\idbm\testData\results\tsv_files\';
        subjId=get(handles.edit1,'String');
        testId=get(handles.edit2,'String');
        testTypeval = get(handles.popupmenu2,'Value');
        SNRval = get(handles.popupmenuSNR,'Value');
        SNR=-1*(2+SNRval*-2);
        fileName=[resultsPath,...
            num2str(SNR),'-',num2str(testTypeval),'.xls'];            
        data1=get(handles.uitable3,'UserData');
        len=size(data1,1);
        data2=cell(1);
        data3=cell(1);
        data2{1}=subjId;
        data3{1}=testId;
        data2=repmat(data2,[len 1]);
        data3=repmat(data3,[len 1]);        
        data4=[data2,data3,data1];
        data4(1,:) = {'SubjectId','TestId','SentenceId','NoWords',...
            'NoRight','PerRight','Words','WordsRight','SentenceType'};
        cell2csv('results.csv', data4); 
    end
    


% --- Executes on button press in pushbuttonClear.
function pushbuttonClear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 set(handles.uitable3,'UserData',[]);
 set(handles.uitable3,'Data',[]);
 
 
 


% --- Executes on button press in pushbuttonResults.
function pushbuttonResults_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Here we retreive the data stores in uitable3 and plot the average number
% of words identified correctly for each condition in a a graph.
% Author: Arun.P.U.
data2=get(handles.uitable3,'UserData');
data3=sortrows(data2,7);
[a,b,c]=unique(data3(:,7),'last');
number = length(a);
switch number
    case 2
        avg(1)=mean(cell2mat(data3(1:b(1),4)));            
        avg(2)=mean(cell2mat(data3(b(1)+1:b(2),4)));
        axes(handles.axes2)
        bar(avg,0.4);
        xlabel(['IDBM','----','Noisy']);
        pr=avg(2)/avg(1)*100;
        set(handles.text2,'String',['% Noisy/IDBM is ',num2str(pr)]);
    case 5
        avg(1)=mean(cell2mat(data3(1:b(1),4)));            
        avg(2)=mean(cell2mat(data3(b(1)+1:b(2),4)));
        avg(3)=mean(cell2mat(data3(b(2)+1:b(3),4)));
        avg(4)=mean(cell2mat(data3(b(3)+1:b(4),4)));
        avg(5)=mean(cell2mat(data3(b(4)+1:b(5),4)));
        axes(handles.axes2)
        bar(avg,0.2);
        xlabel(['Audio','---','AudioVideo','---','IDBM','---','Noisy','---','Video']);          
end
  
