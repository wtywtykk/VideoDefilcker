close all;
clear;

InFile='Rot_YanHe.avi';
OutFile='YanHe_Deflick';
FrameLimit=100000000;
CancelMinFreq=400;
LumXRange=534:1403;
LumYRange=210:915;
LPFFpass=5;
LPFFstop=6;

% InFile='Rot_BadApple.avi';
% OutFile='BadApple_Deflick';
% FrameLimit=100000000;
% CancelMinFreq=400;
% LumXRange=591:1395;
% LumYRange=234:915;
% LPFFpass=2;
% LPFFstop=3;

if 1==0
    ReaderObj=VideoReader(InFile);
    vidWidth=ReaderObj.Width;
    vidHeight=ReaderObj.Height;
    FrameData=readFrame(ReaderObj);
    figure;
    image(FrameData);
    FrameData=FrameData(LumYRange,LumXRange,:);
    figure;
    image(FrameData);
    clear ReaderObj;
else

    if 1==1
        ReaderObj=VideoReader(InFile);
        vidWidth=ReaderObj.Width;
        vidHeight=ReaderObj.Height;
        FrameData=zeros(vidHeight,vidWidth,3);
        k=1;
        while hasFrame(ReaderObj) && (k<FrameLimit)
            FrameData=readFrame(ReaderObj);
            FrameData=FrameData(LumYRange,LumXRange,:);
            Luminance(k)=sum(sum(sum(FrameData)))/vidWidth/vidHeight/3;
            k=k+1;
        end
        clear ReaderObj;
    end
    %Create frequency estimation table
    WindowFunc=hann(length(Luminance))';
    OverSamplingRatio=4096;
    OverSamplingWindowFunc=WindowFunc;
    OverSamplingWindowFunc(length(Luminance)*OverSamplingRatio)=0;
    OverSamplingFFT=fftshift(fft(OverSamplingWindowFunc));
    OverSamplingFFTMag=abs(OverSamplingFFT);
    OverSamplingIndex=-OverSamplingRatio:OverSamplingRatio;
    OverSamplingFreq=OverSamplingIndex/OverSamplingRatio;
    OffsetCoeff=zeros(size(OverSamplingIndex));
    MagCoeff=zeros(size(OverSamplingIndex));
    PhaCoeff=zeros(size(OverSamplingIndex));
    for ii=1:length(OverSamplingIndex)
        IndexC=OverSamplingIndex(ii)+length(OverSamplingFFTMag)/2;
        OffsetCoeff(ii)=(OverSamplingFFTMag(IndexC-OverSamplingRatio)-OverSamplingFFTMag(IndexC+OverSamplingRatio))/OverSamplingFFTMag(IndexC);
        MagCoeff(ii)=OverSamplingFFTMag(length(OverSamplingFFTMag)/2)/OverSamplingFFTMag(IndexC);
        PhaCoeff(ii)=OverSamplingFFT(IndexC);
    end
    PhaCoeff=angle(PhaCoeff);
    FreqPolyFit=polyfit(OffsetCoeff,OverSamplingFreq,7);
    FreqPolyFit(2)=0;
    FreqPolyFit(4)=0;
    FreqPolyFit(6)=0;
    FreqPolyFit(8)=0;
    MagPolyFit=polyfit(OffsetCoeff,MagCoeff,7);
    MagPolyFit(1)=0;
    MagPolyFit(3)=0;
    MagPolyFit(5)=0;
    MagPolyFit(7)=0;
    PhaPolyFit=polyfit(OffsetCoeff,PhaCoeff,7);
    PhaPolyFit(2)=0;
    PhaPolyFit(4)=0;
    PhaPolyFit(6)=0;
    PhaPolyFit(8)=0;

    LuminanceCali=Luminance;
    SinIndex=0:length(LuminanceCali)-1;
    ModSeq=ones(size(SinIndex));

    for ii=1:50
        LuminanceFFT=fft(LuminanceCali.*WindowFunc);
        LuminanceFFTMag=abs(LuminanceFFT);
        LuminancDC=LuminanceFFTMag(1);
        LuminanceFFTMag(1:CancelMinFreq)=0;
        LuminanceFFTMag=LuminanceFFTMag(1:ceil(length(LuminanceFFTMag)/2));
        [Max,MaxIndex]=max(LuminanceFFTMag);
        
        if MaxIndex>CancelMinFreq+1
            LuminanceOffsetCoeff=(LuminanceFFTMag(MaxIndex-1)-LuminanceFFTMag(MaxIndex+1))/LuminanceFFTMag(MaxIndex);
            Offset=polyval(FreqPolyFit,LuminanceOffsetCoeff);
            Mag=LuminanceFFTMag(MaxIndex)*polyval(MagPolyFit,LuminanceOffsetCoeff);
            Pha=angle(LuminanceFFT(MaxIndex))-polyval(PhaPolyFit,LuminanceOffsetCoeff);
            FreqPos=MaxIndex-Offset;
            Freq=2*pi*(FreqPos-1)/length(LuminanceCali);

            ModSeq=ModSeq+2*Mag/LuminancDC*cos(Freq*SinIndex+Pha);
            LuminanceCali=Luminance./ModSeq;
        end
    end

    LuminanceCaliFiltered=LuminanceCali;
    for ii=2:length(LuminanceCaliFiltered)-1
        AverageLum=(LuminanceCaliFiltered(ii-1)+LuminanceCaliFiltered(ii+1))/2;
        DiffLum=LuminanceCaliFiltered(ii-1)-LuminanceCaliFiltered(ii+1);
        if abs(LuminanceCaliFiltered(ii)-AverageLum)>2 && abs(DiffLum)<2.5
            LuminanceCaliFiltered(ii)=AverageLum;
        end
    end

    Fs=60;
    Dpass=0.0057563991496;
    Dstop=0.001;
    flag='scale';
    [N,Wn,BETA,TYPE]=kaiserord([LPFFpass LPFFstop]/(Fs/2), [1 0], [Dstop Dpass]);
    b=fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
    Hd=dfilt.dffir(b);
    LuminanceCaliFiltered(length(LuminanceCali)+Hd.order/2)=0;
    LuminanceCaliFiltered=filter(Hd,LuminanceCaliFiltered);
    LuminanceCaliFiltered=LuminanceCaliFiltered(Hd.order/2+1:length(LuminanceCaliFiltered));

    CaliSeq=LuminanceCaliFiltered./Luminance;

    figure;
    hold on;
    plot(Luminance);
    hold on;
    plot(LuminanceCaliFiltered);
    figure;
    plot(log(abs(fft(Luminance.*WindowFunc))));
    hold on;
    plot(log(abs(fft(LuminanceCaliFiltered.*WindowFunc))));

    if 1==1
        ReaderObj=VideoReader(InFile);
        WriterObj=VideoWriter(OutFile,'MPEG-4');
        WriterObj.Quality=100;
        WriterObj.FrameRate=ReaderObj.FrameRate;
        open(WriterObj);
        vidWidth=ReaderObj.Width;
        vidHeight=ReaderObj.Height;
        FrameData=zeros(vidHeight,vidWidth,3);
        k=1;
        while hasFrame(ReaderObj) && (k<FrameLimit)
            FrameData=readFrame(ReaderObj);
            FrameData=FrameData*CaliSeq(k);
            writeVideo(WriterObj,FrameData);
            100*k/length(CaliSeq)
            k=k+1;
        end
        close(WriterObj);
    end

    clear ReaderObj;
    clear WriterObj;

end
