program NewPATC052115BalancedLongInt;

{Note this is a very extensive program that mostly does one calculation for a sliding window of sequence then keeps track of overall scores while
the window slides.  The key parts of the algorithm are contained between lines 1442 and 1575, the rest is all fancy book-keeping}
{Important variables}
   Const {Untyped Constants}
     InputFileBase='patc.bat';
     Version='PATCV';
     MoveNumber=3; {Number of moves ahead the program will look for optimizing bend element}
     FourToTheMoveNumber=64; {4^MoveNumber}
     FiveBondProfit=30;
     FourBondProfit=20;
     ThreeBondProfit=10;
     TwoBondProfit=0;
     NoBondProfit=-5;
     NineBaseLoss=16;
     TenBaseLoss=8;
     ElevenBaseLoss=16;
     TwelveBaseLoss=32;
     ActiveWindow=1200;
     ActiveWindowMinusOne=1199;
     TwoTimesActiveWindow=1280;
     TwoTimesActiveWindowMinusOne=1279;
     HistogramMax=1024;
     CPMax=1024;
     BestBendsRetained=10;
     Pi=3.1415926535;

   Type
     Possibility = Record
        MArray:Array[1..MoveNumber] of integer;
        Cost:Array[1..MoveNumber] of integer;
     end;

   CharArray = Array [0..TwoTimesActiveWindowMinusOne] of char;
   BoolArray = Array [0..TwoTimesActiveWindowMinusOne] of Boolean;
   IntegerArray = Array [0 .. TwoTimesActiveWindowMinusOne] of Integer;
   LongIntegerArray = Array [0 .. TwoTimesActiveWindowMinusOne] of Int64;
   PossibilityArrayT = Array [1.. FourToTheMoveNumber] of Possibility;
   FiveBaseArray = Array [0 .. 1023] of integer;
   FiveBaseLongArray = Array [0 .. 1023] of Int64;
   ShortPossibilityArrayT = Array [1..20] of integer;
   StringArray = Array [1..BestBendsRetained] of String;
   BestArray = Array [1..BestBendsRetained] of Int64;
   HistogramArray = Array [0..HistogramMax] of Int64;
   CPArray = Array [0..CPMax] of Int64;
   PArray = Array [0..14,0..4] of Int64;
   PArray2 = Array [0..14,0..16] of Int64;

Var
    RecordSeparationArray, RecordSNPArray, ThoroughSeparationArray, CenterOff, SeparationArrayForLowerOnly, SeparationArrayForExonJumpover, SeparationArrayNoExonJumpover, IsUpper, FoundUpper, RecordOn: Boolean;
    {SumaryStats}
    AllfN,AllfP,AlluN,AlluP,AlleN,AlleP,AlliN,AlliP,AlldN,AlldP,AllhN,AllhP,lowertrail: Int64;
    SuperSensitive: Boolean; {This tells the program to start looking at single AA/TT basepairs as significant}
    Nxde:Integer;
    Typearray:CharArray;
    CurType,OldType:Char;
    CurrentRepresentation: LongInt; {The current representation; *4, mod 1024 + next base to get next representation}
        {The meaning of a representation in base sequence is a five base sequence in 5'->3' orientation stopping at the current base}
    Iade:longint;
    Xbase1, Xbase2 ,AATTCount,ATTACount,GCount, ACount, TCount, CCount, NCount,XCount,Dashcount,Curpair,ExternalIndex,SubSequenceStart,Direction,SteinLength,FeatureLength:Int64;
    {Curpair is a 2-base segment equivalent to current representation}
    RandomCharacter, Scrambler: Integer;
    Threshold,IntronF,UpsF,DwnF,ExonF,FiveF,ThreeF,IntronAF,ExonAF,FiveAF,ThreeAF,IntronD,UpsD,DwnD,ExonD,FiveD,ThreeD,IntronAD,ExonAD,FiveAD,ThreeAD:Double;
    Number1,IntronP,UpsP,DwnP,ExonP,FiveP,ThreeP,IntronAP,ExonAP,FiveAP,ThreeAP: Double;
    IntronB,UpsB,DwnB,ExonB,FiveB,ThreeB,IntronAB,ExonAB,FiveAB,ThreeAB,Goodbase: Int64;
    IntronN,UpsN,DwnN,ExonN,FiveN,ThreeN,IntronAN,ExonAN,FiveAN,ThreeAN: Int64;
    UpsQ,DwnQ,IntronQ,ExonQ,FiveQ,ThreeQ,IntronAQ,ExonAQ,FiveAQ,ThreeAQ,SteinFormat: Boolean;
    FiveFiveAN,FiveGeneN,GeneAGeneN,GeneThreeN,ThreeAThreeN: Int64;
    FileNameBase,SequenceName,SubSequenceName,BinType,BaseString,Chromosome,CommonName,WorkingDirectory:String;
        InputFileName,SeqFileName,SeqFileName0,SeqFileFormat,OFileName,DataFileName,DisFileName,NumFileName, SourceFileName,Curline,LSS,Description,TheDate,TheTime,ParamaterList:String;
        InputFile,SeqFile,DataFile,OFile,DisFile,SourceFile,NumFile:Text;
        Tab: Char;
    RepresentationArray,AddArray: FiveBaseArray; {Benefits of each representation}

    PhaseArray: PArray; {Base composition at each point in the phasing cycle, 14 is for unphased DNA}
    PhaseArray2: PArray2; {Base composition at each point in the phasing cycle, 14 is for unphased DNA}
    OffBaseArray: FiveBaseLongArray; {sequences on the opposite face of the helix}
    OnBaseArray: FiveBaseLongArray; {sequences on this face of the helix}

    ShortPossibilityArray:ShortPossibilityArrayT; {Summary of possibility costs}

    UxA,FxA,IxA,ExA,HxA,DxA,UnA,FnA,InA,EnA,HnA,DnA,InA50,IxA50,EnA50,ExA50:HistogramArray; {Histograms based on positions within sequece elements [stein format only]}
    Tn1,Tn2,Tn3,Tn4:Int64; {Tn1 describes the distance from the start/end of the gene, Tn2 runs a 1-50 loop, and Tn3 is Tr1 div 50}
    TempRatio:Double;
    Histogram,LookBackArray:HistogramArray; {Global Histogram of values}
    RetroChar,OldRetroChar:Char;
    SequenceTextArray: CharArray; {Window length base sequence}
    SequenceTypeArray: BoolArray; {true for bone fide lower case bases}
    SequenceArray: IntegerArray; {Numerical Sequence values}
    RArray: IntegerArray; {Representations for each base}
    ValueArray,ValueArrayB,SepValueArray: IntegerArray; {Values for Inflexibility in 5-base interval centered on that base}
    Choicearray: IntegerArray; {Choices for backdraft that have already been calculated}
    BestFitArray: IntegerArray; {Keeps track of best fit value for feature}
    BestFitPhase: IntegerArray; {Keeps track of the phasing of this base relative to the best fit for that position, Phase=0 means this is the peak of phasing}
    PhasedN: Int64; {This is the number of phased position in the last window of Analysis}
    PossibilityArray: PossibilityArrayT; {Array of possibilities}
    SeparationArray,SeparationNArray,SNPArray: LongIntegerArray; {Distribution of Separations between phasing peaks}
    CurPhase: Int64; {The current phase in the helix: 0 is the end of an A-track}
    CurPeriod: Int64; {The current periodicity of the helix: should be 9-12}
    PosG: Int64; {Base Number in Genome}
    PosC: Int64; {Base Number in Chromosome}
    PosW: Int64; {Base number in whole sequence}
    PosS: Int64; {Base Number in Subsequence}
    PosT: Int64; {Number of bases tested in current stretch}
    PosR: Int64; {Number of bases recorded in bin}
    PosB: Int64; {Number of bases recorded in current stretch}
    PosL: Int64; {Base Number within Line}
    PosA: Int64; {Evaluate Position within 512b Cycle}
    PosD: Int64; {Record Position within 512b Cycle}
    PosX: Int64; {Position within a worm autosome}
    cub,cup,cuu,nub,nup,nuu,rb,adb: Int64;
    {adb=all determined bases in bin}
    {rb= all "x" bases in bin}
    {cub= all capitalized bases}
    {cup= phased capitalized bases}
    {cuu= unphased capitailized bases}
    {nub= all non-capitalized bases}
    {nup= phased non-capitalized bases}
    {nuu= unphased non-capitailized bases}
    PosNR,PosNA: Int64; {Index for separation list NR: REleative, NA: Absolute}
    SeqFileNum: Int64; {Autosome we're currently looking at}
    BinSize: Int64; {Bin Size}
    BinNumber: Int64; {Bin Number}
    SubSequenceNumber: Int64; {Bin Number}
    BendsInBin: Int64; {Bends in this bin}
    PointsInBin:Double; {Total Bend Points in this bin}
    BendFrequency: Double; {Bendsinbin/Binsize}
    BendDensity: Double; {Pointsinbin/Binsize}
    TotalBends,SubBends: Int64;
    TotalPoints,SubPoints: Double;
    TotalBases,SubBases: Int64; {Total bases examined}
    TotalBendFrequency,SubBendFrequency: Double;
    TotalPointDensity,SubPointDensity: Double;
    UserStillGoing,AutoRun: Boolean;
    BestScores,BestPositions:BestArray;
    BestStrings:StringArray;
    WorstBestScore,WorstBest:Int64;
    ln2:double;
    CPA,CTA,TPA,TTA:CPArray;
    nTA:Int64;
    ReflectCPArray:Boolean;
    {CPA is a histogram of bend incidences as a function of distance from the end of the current subsequence}
    {CTA is a histogram of total bases as a function of distances from the end of the current subsequence}
    {TPA is a compendium of bend incidences as a function of distance from the end of the current subsequence}
    {CTA is a compendium of total bases as a function of distances from the end of the current subsequence}
    {nTA is the number of bases n each bin}
    NumericalArray:Boolean;
    {Shall we make an extra output file (*.num) with structure:
       > Gene Name
       Base Number <tab> Phase [0-13] <tab> Phasing Score
       > Next Gene Name
       ...
       This is for graphic output of individual genes}
    FilterArray:Boolean; {A second use for the .num file, in this case a filtered
        set of sequences consisting of all bases for which the cutoff value is bested.  Other bases are all
        N's in this sequence}
    FilterCount:Integer; {How many characters in the current output line}
    fws1,fwe1,fwc1,fwz1:Int64; {start, end, and current indexes relative to the filter termini}
    fwb1:boolean; {should we record this base}
    BalanceMode:boolean; {?subtract off the any contributions for all contributing AA/TT dinucleotides on the wrong side of helix?
    						Hoping this will reduce spurious PATC calls in highly AT-rich genomes and genomic regions}

Procedure CPDouble;  var icpe:longint;
begin
   nTA:=nTA*2;
   icpe:=0;
   repeat
     CPA[icpe]:=CPA[2*icpe]+CPA[2*icpe+1];
     CTA[icpe]:=CTA[2*icpe]+CTA[2*icpe+1];
     TPA[icpe]:=TPA[2*icpe]+TPA[2*icpe+1];
     TTA[icpe]:=TTA[2*icpe]+TTA[2*icpe+1];
     icpe:=icpe+1;
   until icpe=CPMax div 2;
   repeat
     CPA[icpe]:=0;
     CTA[icpe]:=0;
     TPA[icpe]:=0;
     TTA[icpe]:=0;
     icpe:=icpe+1;
   until icpe=CPMax+1;
end;

Procedure CPClear;  var icpe:longint;
begin
   for icpe:=0 to CPMax do begin
     CPA[icpe]:=0;
     CTA[icpe]:=0;
   end;
end;

Procedure TPClear;  var icpe:longint;
begin
   nTA:=1;
   for icpe:=0 to CPMax do begin
     TPA[icpe]:=0;
     TTA[icpe]:=0;
   end;
end;

Procedure CPRecord;  var icpe,jcpe:longint;
begin
   for icpe:=0 to CPMax do begin
     jcpe:=icpe;
     if icpe*nTA*2>PosS then begin
        If ReflectCPArray then jcpe:=(SubBases-(icpe*nTA)) div nTA;
        If Not(ReflectCPArray) then jcpe:=icpe;
     end; {if icpe*nTA=}
     if jcpe>-1 then begin
        TPA[jcpe]:=TPA[jcpe]+CPA[icpe];
        TTA[jcpe]:=TTA[jcpe]+CTA[icpe];
     end; {if jcpe>-1}
   end;
end;

Procedure BinEndHandle;
begin
   BendFrequency:=0;
   BendDensity:=0;
  If ACount+TCount+GCount+CCount>0 then begin
      BendFrequency:=BendsInBin/(ACount+TCount+GCount+CCount);
      BendDensity:=PointsInBin/(ACount+TCount+GCount+CCount);
   end;
   IntronF:=0; IntronD:=0; If IntronN>0 then begin IntronF:=IntronB/IntronN; IntronD:=IntronP/IntronN; end;
   FiveF:=0; FiveD:=0; If FiveN>0 then begin FiveF:=FiveB/FiveN; FiveD:=FiveP/FiveN; end;
   ThreeF:=0; ThreeD:=0; If ThreeN>0 then begin ThreeF:=ThreeB/ThreeN; ThreeD:=ThreeP/ThreeN; end;
   ExonF:=0; ExonD:=0; If ExonN>0 then begin ExonF:=ExonB/ExonN; ExonD:=ExonP/ExonN; end;
   upsF:=0; upsD:=0; If upsN>0 then begin upsF:=upsB/upsN; upsD:=upsP/upsN; end;
   dwnF:=0; dwnD:=0; If dwnN>0 then begin dwnF:=dwnB/dwnN; dwnD:=dwnP/dwnN; end;

   IntronAF:=0; IntronAD:=0; If IntronAN>0 then begin IntronAF:=IntronAB/IntronAN; IntronAD:=IntronAP/IntronAN; end;
   FiveAF:=0; FiveAD:=0; If FiveAN>0 then begin FiveAF:=FiveAB/FiveAN; FiveAD:=FiveAP/FiveAN; end;
   ThreeAF:=0; ThreeAD:=0; If ThreeAN>0 then begin ThreeAF:=ThreeAB/ThreeAN; ThreeAD:=ThreeAP/ThreeAN; end;
   ExonAF:=0; ExonAD:=0; If ExonAN>0 then begin ExonAF:=ExonAB/ExonAN; ExonAD:=ExonAP/ExonAN; end;
  If CommonName<>'NoCommonName' then Write(Datafile,CommonName) else Write(Datafile,SubSequenceName);
   Write(Datafile,Tab);
   Write(Datafile,SubSequenceName);
   Write(Datafile,Tab);
   Write(Datafile,IntronF);
   Write(Datafile,Tab);
   Write(Datafile,ExonF);
   Write(Datafile,Tab);
   Write(Datafile,FiveF);
   Write(Datafile,Tab);
   Write(Datafile,ThreeF);
   Write(Datafile,Tab);
   Write(Datafile,IntronD);
   Write(Datafile,Tab);
   Write(Datafile,ExonD);
   Write(Datafile,Tab);
   Write(Datafile,FiveD);
   Write(Datafile,Tab);
   Write(Datafile,ThreeD);
   Write(Datafile,Tab);
   Write(Datafile,IntronB);
   Write(Datafile,Tab);
   Write(Datafile,ExonB);
   Write(Datafile,Tab);
   Write(Datafile,FiveB);
   Write(Datafile,Tab);
   Write(Datafile,ThreeB);
   Write(Datafile,Tab);
   Write(Datafile,IntronN);
   Write(Datafile,Tab);
   Write(Datafile,ExonN);
   Write(Datafile,Tab);
   Write(Datafile,FiveN);
   Write(Datafile,Tab);
   Write(Datafile,ThreeN);
   Write(Datafile,Tab);
   Write(Datafile,IntronP);
   Write(Datafile,Tab);
   Write(Datafile,ExonP);
   Write(Datafile,Tab);
   Write(Datafile,FiveP);
   Write(Datafile,Tab);
   Write(Datafile,ThreeP);
   Write(Datafile,Tab);
   Write(Datafile,IntronAF);
   Write(Datafile,Tab);
   Write(Datafile,ExonAF);
   Write(Datafile,Tab);
   Write(Datafile,FiveAF);
   Write(Datafile,Tab);
   Write(Datafile,ThreeAF);
   Write(Datafile,Tab);
   Write(Datafile,IntronAD);
   Write(Datafile,Tab);
   Write(Datafile,ExonAD);
   Write(Datafile,Tab);
   Write(Datafile,FiveAD);
   Write(Datafile,Tab);
   Write(Datafile,ThreeAD);
   Write(Datafile,Tab);
   Write(Datafile,IntronAB);
   Write(Datafile,Tab);
   Write(Datafile,ExonAB);
   Write(Datafile,Tab);
   Write(Datafile,FiveAB);
   Write(Datafile,Tab);
   Write(Datafile,ThreeAB);
   Write(Datafile,Tab);
   Write(Datafile,IntronAN);
   Write(Datafile,Tab);
   Write(Datafile,ExonAN);
   Write(Datafile,Tab);
   Write(Datafile,FiveAN);
   Write(Datafile,Tab);
   Write(Datafile,ThreeAN);
   Write(Datafile,Tab);
   Write(Datafile,IntronAP);
   Write(Datafile,Tab);
   Write(Datafile,ExonAP);
   Write(Datafile,Tab);
   Write(Datafile,FiveAP);
   Write(Datafile,Tab);
   Write(Datafile,ThreeAP);
   Write(Datafile,Tab);
   Write(Datafile,SubSequenceNumber);
   Write(Datafile,Tab);
   Write(Datafile,BinNumber);
   Write(Datafile,Tab);
   Write(Datafile,BinType);
   Write(Datafile,Tab);
   Write(Datafile,PosR);
   Write(Datafile,Tab);
   Write(Datafile,PosS-PosR+1);
   Write(Datafile,Tab);
   Write(Datafile,PosW-Direction*PosR);
   Write(Datafile,Tab);
   Write(Datafile,GCount);
   Write(Datafile,Tab);
   Write(Datafile,ACount);
   Write(Datafile,Tab);
   Write(Datafile,TCount);
   Write(Datafile,Tab);
   Write(Datafile,CCount);
   Write(Datafile,Tab);
   Write(Datafile,NCount);
   Write(Datafile,Tab);
   Write(Datafile,XCount);
   Write(Datafile,Tab);
   Write(Datafile,Dashcount);
   Write(Datafile,Tab);
   Write(Datafile,AATTCount);
   Write(Datafile,Tab);
   Write(Datafile,ATTACount);
   Write(Datafile,Tab);
   Write(Datafile,Direction);
   Write(Datafile,Tab);
   Write(Datafile,BendsInBin);
   Write(Datafile,Tab);
   Write(Datafile,BendFrequency);
   Write(Datafile,Tab);
   Write(Datafile,PointsInBin);
   Write(Datafile,Tab);
   Write(Datafile,BendDensity);
   Write(Datafile,Tab);
   Write(Datafile,cub);
   Write(Datafile,Tab);
   Write(Datafile,cup);
   Write(Datafile,Tab);
   Write(Datafile,cuu);
   Write(Datafile,Tab);
   Write(Datafile,nub);
   Write(Datafile,Tab);
   Write(Datafile,nup);
   Write(Datafile,Tab);
   Write(Datafile,nuu);
   Write(Datafile,Tab);
   Write(Datafile,rb);
   Write(Datafile,Tab);
   Write(Datafile,adb);
   Write(Datafile,Tab);
      Write(Datafile,upsF);
   Write(Datafile,Tab);
   Write(Datafile,dwnF);
   Write(Datafile,Tab);
    Write(Datafile,exp(upsD*ln2/10));
   Write(Datafile,Tab);
   Write(Datafile,exp(dwnD*ln2/10));
   Write(Datafile,Tab);
   Write(Datafile,upsB);
   Write(Datafile,Tab);
   Write(Datafile,dwnB);
   Write(Datafile,Tab);
   Write(Datafile,upsN);
   Write(Datafile,Tab);
   Write(Datafile,dwnN);
   Write(Datafile,Tab);
   Write(Datafile,upsD);
   Write(Datafile,Tab);
   Write(Datafile,dwnD);
   Write(Datafile,Tab);
   Write(Datafile,ExternalIndex);
   Write(Datafile,Tab);
   Writeln(Datafile);
   Write(PosS);
   Write(Tab);
   Write(adb);
   Write(Tab);
   Write(BendsInBin);
   Write(Tab);
   Write(BendFrequency);
   Writeln;
   BendsInBin:=0;
   PointsInBin:=0;
   UpsB:=0;DwnB:=0;IntronB:=0;ExonB:=0;ThreeB:=0;FiveB:=0;
   UpsN:=0;DwnN:=0;IntronN:=0;ExonN:=0;ThreeN:=0;FiveN:=0;
   UpsP:=0;DwnP:=0;IntronP:=0;ExonP:=0;ThreeP:=0;FiveP:=0;
   IntronAB:=0;ExonAB:=0;ThreeAB:=0;FiveAB:=0;
   IntronAN:=0;ExonAN:=0;ThreeAN:=0;FiveAN:=0;
   IntronAP:=0;ExonAP:=0;ThreeAP:=0;FiveAP:=0;
   PosR:=0;
   BinNumber:=BinNumber+1;
   AATTCount:=0;
   ATTACount:=0;
   GCount:=0;
   ACount:=0;
   TCount:=0;
   CCount:=0;
   NCount:=0;
   XCount:=0;
   Dashcount:=0;
   cub:=0;
   cup:=0;
   cuu:=0;
   nub:=0;
   nup:=0;
   nuu:=0;
   adb:=0;
   rb:=0;
end; {BinEndHandle}

Procedure RecordStatus;
var
   cvde,bffe,hffe,nffe,fwi1:longint;
   rffe:double;
begin
         RetroChar:=SequenceTextArray[PosD];
         CurPair:=RArray[PosD] mod 16;
         hffe:=BestFitArray[PosD];
         If hffe>HistogramMax then hffe:=HistogramMax;
         If hffe<0 then hffe:=0;
         Histogram[hffe]:=Histogram[hffe]+1;
         PointsInBin:=PointsInBin+BestFitarray[PosD];
         TotalPoints:=TotalPoints+BestFitarray[PosD];
         SubPoints:=SubPoints+BestFitarray[PosD];
         bffe:=14; {assume first no phasing}
         If PosS>nTA*CPMax then CPdouble;
         if nTA>0 then CTA[PosS div nTA]:=CTA[PosS div nTA]+1;

        { A debug routine:
        Not needed for standard program use
        If PosD=221 then begin
              For cvde:=0 to 1280 do begin
             Write(SequenceTextArray[cvde]);
             If cvde mod 50=0 then Writeln;
           end;
           ReadLn;
         end;}

         If BestFitarray[PosD]>Threshold then begin
            BendsInBin:=BendsInBin+1; TotalBends:=TotalBends+1; SubBends:=SubBends+1;
            PhasedN:=PhasedN+1;
            if nTA>0 then CPA[PosS div nTA]:=CPA[PosS div nTA]+1;
            bffe:=BestFitPhase[PosD];
            If bffe=0 then OnBaseArray[RArray[PosD]]:=OnBaseArray[RArray[PosD]]+1;
            If bffe=5 then OffBaseArray[RArray[PosD]]:=OffBaseArray[RArray[PosD]]+1;
            If bffe<0 then bffe:=0;
            If bffe>13 then bffe:=13;

            end; {If BestFitarray[PosD]>Threshold}
            If RecordSNPArray AND (PosD<512) AND (PosD>0) then begin
            SNPArray[PosD]:=SNPArray[PosD]+ValueArray[PosD];
            {If Retrochar='M' then writeln(PosD);}
            If PosD=250 then begin
            	Nxde:=601;
           		Case RetroChar of
            		'R','r','Y','y':Nxde:=700;
            		'K','k','M','m':Nxde:=750;
            		'A','a','T','t':Nxde:=800;
            		'G','g','C','c':Nxde:=850;
            		'S','s':Nxde:=900;
            		'W','w':Nxde:=950;
            	end; {case}
            If (BestFitPhase[PosD]>-1) and (BestFitPhase[PosD]<15) then Nxde:=Nxde+BestFitPhase[PosD]+1;
            SNPArray[Nxde]:=SNPArray[Nxde]+1;
            end; {250}

            end; {RecordSNPArray}
            If NumericalArray then begin
               Write(Numfile,PosS+1);
               Write(Numfile,tab);
               Write(Numfile,RetroChar);
               Write(Numfile,tab);
               Write(Numfile,BestFitArray[PosD]);
               Write(Numfile,tab);
               Write(Numfile,BestFitPhase[PosD]);
               Write(Numfile,tab);
               nffe:=BestFitPhase[PosD];
               rffe:=Sin(nffe*2*pi/10.2);
               If nffe=-999 then rffe:=0;
               Write(Numfile,rffe);
               Writeln(Numfile,tab);
            end;
            If FilterArray then begin
               fwb1:=false;
               fwc1:=fwc1+1;
               If BestFitArray[PosD]>threshold then fwc1:=0;
               If fwc1<fws1+1 then fwb1:=true;
               If fwc1>fws1+100 then fwc1:=fws1+100;
               fwz1:=PosD;
               For fwi1:=1 to fwe1 do begin
                  fwz1:=fwz1+1;
                  if fwz1=TwoTimesActiveWindow then fwz1:=0;
                  If BestFitArray[fwz1]>threshold then fwb1:=true;
               end; {for fwi1=}
               If fwb1 then Write(Numfile,RetroChar) else Write(Numfile,'n');
               FilterCount:=FilterCount+1;
               If FilterCount>60  Then begin Writeln(Numfile); FilterCount:=1 end;
            end;

         PosW:=PosW+Direction;
         PosR:=PosR+1;
         PosB:=PosB+1;
         PosS:=PosS+1;
         Case RetroChar of
           'G','g': begin GCount:=GCount+1; TotalBases:=TotalBases+1; SubBases:=SubBases+1; PhaseArray[bffe,0]:=PhaseArray[bffe,0]+1; PhaseArray2[bffe,CurPair]:=PhaseArray2[bffe,CurPair]+1;PhaseArray2[bffe,16]:=PhaseArray2[bffe,16]+1;end;
           'A','a': begin ACount:=ACount+1; TotalBases:=TotalBases+1; SubBases:=SubBases+1; PhaseArray[bffe,1]:=PhaseArray[bffe,1]+1; PhaseArray2[bffe,CurPair]:=PhaseArray2[bffe,CurPair]+1;PhaseArray2[bffe,16]:=PhaseArray2[bffe,16]+1;end;
           'T','t': begin TCount:=TCount+1; TotalBases:=TotalBases+1; SubBases:=SubBases+1; PhaseArray[bffe,2]:=PhaseArray[bffe,2]+1; PhaseArray2[bffe,CurPair]:=PhaseArray2[bffe,CurPair]+1;PhaseArray2[bffe,16]:=PhaseArray2[bffe,16]+1;end;
           'C','c': begin CCount:=CCount+1; TotalBases:=TotalBases+1; SubBases:=SubBases+1; PhaseArray[bffe,3]:=PhaseArray[bffe,3]+1; PhaseArray2[bffe,CurPair]:=PhaseArray2[bffe,CurPair]+1;PhaseArray2[bffe,16]:=PhaseArray2[bffe,16]+1;end;
           'N','n': begin NCount:=NCount+1; PhaseArray[bffe,4]:=PhaseArray[bffe,4]+1; end;
           'X','x': begin XCount:=XCount+1; PhaseArray[bffe,4]:=PhaseArray[bffe,4]+1; end;
           '-': begin DashCount:=DashCount+1; PhaseArray[bffe,4]:=PhaseArray[bffe,4]+1; end;
         end;
         Case RetroChar of
           'G','A','T','C': begin
              cub:=cub+1;
              adb:=adb+1;
              If BestFitarray[PosD]>Threshold then cup:=cup+1 else cuu:=cuu+1;
              TypeArray[PosD]:='e';
              If OldType='f' then begin Tn1:=0;Tn2:=0;Tn3:=0; end;
              ExonN:=ExonN+1;
              ExonP:=ExonP+BestFitarray[PosD];
              if SteinLength=0 then SteinLength:=1;
              Tn4:=((Tn1*200) div SteinLength)+100;
             { Debug routine Write('Tn1:');Writeln(Tn1);  Write('SteinLength:');Writeln(SteinLength);}
              If Tn4>300 then Tn4:=300;
              If BestFitarray[PosD]>Threshold then begin
               ExonB:=ExonB+1;
               AlleP:=AlleP+1;
               ExA[Tn4]:=ExA[Tn4]+1;
               ExA50[Tn3]:=ExA50[Tn3]+1;
              end; {>Threshold}
               EnA[Tn4]:=EnA[Tn4]+1;
               AlleN:=AlleN+1;
               EnA50[Tn3]:=EnA50[Tn3]+1;
               Tn2:=Tn2+1; If Tn2=50 then begin Tn2:=0; If Tn3<HistogramMax then Tn3:=Tn3+1; end;
               Tn1:=Tn1+1;
             If ExonN<201 then begin
             	 ExonAN:=ExonAN+1;
            	 ExonAP:=ExonAP+BestFitarray[PosD];
           		 If BestFitarray[PosD]>Threshold then ExonAB:=ExonAB+1;
           	  end;
              end; {GATC}
           'g','a','t','c': begin
               nub:=nub+1;
               adb:=adb+1;
               If BestFitarray[PosD]>Threshold then nup:=nup+1 else nuu:=nuu+1;
               If ExonN=0 then begin
                 If TypeArray[PosD]='f' then begin
                    If OldType<>'f' then begin Tn1:=0;Tn2:=0;Tn3:=0; end;
                    FiveN:=FiveN+1;
                    FiveP:=FiveP+BestFitarray[PosD];
                    Tn4:=Tn1 div 10;
                    If Tn4>99 then Tn4:=99;
                    If BestFitarray[PosD]>Threshold then begin
                        FxA[Tn4]:=FxA[Tn4]+1; FiveB:=FiveB+1;
						If Tn1>300 then AllfP:=AllfP+1;
                    end; {>Threshold}
				  If Tn1>300 then AllfN:=AllfN+1;
                  FnA[Tn4]:=FnA[Tn4]+1;
                  Tn1:=Tn1+1;
                    If FiveN>800 then begin
             	       FiveAN:=FiveAN+1;
            	       FiveAP:=FiveAP+BestFitarray[PosD];
           		       If BestFitarray[PosD]>Threshold then FiveAB:=FiveAB+1;
           	        end;
           	     end; {f}
                 If TypeArray[PosD]='u' then begin
                 If OldType='f' then begin Tn1:=0;Tn2:=0;Tn3:=0; end;
                    upsN:=upsN+1;
                    upsP:=upsP+BestFitarray[PosD];
                   if SteinLength=0 then SteinLength:=1;
                   Tn4:=((Tn1*200) div SteinLength)+100;
                   If Tn4>300 then Tn4:=300;
              If BestFitarray[PosD]>Threshold then begin
                   upsB:=upsB+1;
                   If Tn1<200 then begin UxA[Tn4]:=UxA[Tn4]+1; AlluP:=AlluP+1; end;
               end; {>Threshold}
                 If Tn1<200 then begin UnA[Tn4]:=UnA[Tn4]+1; AlluN:=AlluN+1; end;
                 Tn1:=Tn1+1;
                 Tn2:=Tn2+1; If Tn2=50 then begin Tn2:=0; If Tn3<HistogramMax then Tn3:=Tn3+1; end;

           	     end; {f}
           	   end; {ExonN=0}
           	   If Not(SteinFormat) then begin
                If (ExonN>0) AND (PosS<GeneThreeN+1) then begin
                 IntronN:=IntronN+1;
                 IntronP:=IntronP+BestFitarray[PosD];
                 If BestFitarray[PosD]>Threshold then IntronB:=IntronB+1;
                 If IntronN<200 then begin
             	    IntronAN:=IntronAN+1;
            	    IntronAP:=IntronAP+BestFitarray[PosD];
           		    If BestFitarray[PosD]>Threshold then IntronAB:=IntronAB+1;
           	     end; {intronN<200}
           	   end; {ExonN>0 and}
              If (ExonN>0) AND (PosS>GeneThreeN) then begin
                 ThreeN:=ThreeN+1;
                 ThreeP:=ThreeP+BestFitarray[PosD];
                 If BestFitarray[PosD]>Threshold then ThreeB:=ThreeB+1;
                 If ThreeN<200 then begin
             	    ThreeAN:=ThreeAN+1;
            	    ThreeAP:=ThreeAP+BestFitarray[PosD];
           		    If BestFitarray[PosD]>Threshold then ThreeAB:=ThreeAB+1;
           	     end; {IfThreen<}
           	   end; {If ExonN>0}
           	   end {If not Steinformat}
           	   else {if steinformat} begin
           	     if TypeArray[PosD]='i' then begin
           	         IntronN:=IntronN+1;
                     IntronP:=IntronP+BestFitarray[PosD];
                    if SteinLength=0 then SteinLength:=1;
                    Tn4:=((Tn1*200) div SteinLength)+100;
              		If Tn4>300 then Tn4:=300;
					If BestFitarray[PosD]>Threshold then begin
                         IntronB:=IntronB+1;
                          AlliP:=AlliP+1;
                         IxA[Tn4]:=IxA[Tn4]+1;
                         IxA50[Tn3]:=IxA50[Tn3]+1;
                         end;
                          AlliN:=AlliN+1;
                         InA[Tn4]:=InA[Tn4]+1;
                         InA50[Tn3]:=InA50[Tn3]+1;
                         Tn1:=Tn1+1;
                         Tn2:=Tn2+1; If Tn2=50 then begin Tn2:=0; If Tn3<HistogramMax then Tn3:=Tn3+1; end;
                    If IntronN<200 then begin
             	       IntronAN:=IntronAN+1;
            	       IntronAP:=IntronAP+BestFitarray[PosD];
           		       If BestFitarray[PosD]>Threshold then IntronAB:=IntronAB+1;
           	        end;
           	      end;
           	     if TypeArray[PosD]='h' then begin
           	         If OldType<>'h' then begin Tn1:=0;Tn2:=0;Tn3:=0; end;
           	         ThreeN:=ThreeN+1;
                     ThreeP:=ThreeP+BestFitarray[PosD];
                      Tn4:=(Tn1 div 10) + 301;
                      If Tn4>400 then Tn4:=400;
                      If BestFitarray[PosD]>Threshold then begin
                        If Tn1<700 then AllhP:=AllhP+1;
                       HxA[Tn4]:=HxA[Tn4]+1;
                        ThreeB:=ThreeB+1;
                    end; {>Threshold}
                  HnA[Tn4]:=HnA[Tn4]+1;
                  If Tn1<700 then AllhN:=AllhN+1;
                  Tn1:=Tn1+1;
                     If ThreeN<200 then begin
             	       ThreeAN:=ThreeAN+1;
            	       ThreeAP:=ThreeAP+BestFitarray[PosD];
           		       If BestFitarray[PosD]>Threshold then ThreeAB:=ThreeAB+1;
           	        end;
           	      end;
           	     if TypeArray[PosD]='d' then begin
           	         dwnN:=dwnN+1;
                     dwnP:=dwnP+BestFitarray[PosD];
                     if SteinLength=0 then SteinLength:=1;
                     Tn4:=((Tn1*200) div SteinLength)+100;
            	     If Tn4>300 then Tn4:=300;
					 If BestFitarray[PosD]>Threshold then begin
                        dwnB:=dwnB+1;
                        If Tn1>SteinLength-400 then begin DxA[Tn4]:=DxA[Tn4]+1;AlldP:=AlldP+1; end
                      end;
                         If Tn1>SteinLength-400 then begin AlldN:=AlldN+1; DnA[Tn4]:=DnA[Tn4]+1; end;
                          Tn1:=Tn1+1;
                          Tn2:=Tn2+1; If Tn2=50 then begin Tn2:=0; If Tn3<HistogramMax then Tn3:=Tn3+1; end;
           	      end;
               end; {if not stein}
             end;  {gatc}
           'X','x': begin
              adb:=adb+1; rb:=rb+1;
              end;  {Xx}
           'N','n','-': begin end;
         end;
         If RetroChar='a' then RetroChar:='A';
         If RetroChar='t' then RetroChar:='T';
         If ((RetroChar='A') and (OldRetroChar='A')) or ((RetroChar='T') and (OldRetroChar='T')) then AATTCount:=AATTCount+1;
         If ((RetroChar='A') and (OldRetroChar='T')) or ((RetroChar='T') and (OldRetroChar='A')) then ATTACount:=ATTACount+1;
         OldRetroChar:=RetroChar;
         OldType:=TypeArray[PosD];
         If PosR=Binsize then BinEndHandle;
end;

Procedure PurgeArrays;
begin
         If PosT<ActiveWindow then begin
         PosD:=-1;
         While PosB<PosT do begin
            PosD:=PosD+1;
            RecordStatus;
         end; {While PosB<PosT}
         end else
         PosD:=PosA-ActiveWindow;
         If PosD<0 then PosD:=PosD+TwoTimesActiveWindow;
         While PosB<PosT do begin
            PosD:=PosD+1;
            If PosD=TwoTimesActiveWindow then PosD:=0;
            RecordStatus;
         end; {PosT<ActiveWindow}
         If PosR>0 then BinEndHandle;
end;

Procedure BinClear;
var inbe:longint;
begin
     For Inbe:=0 to TwoTimesActiveWindowMinusOne do begin
             SequenceTextArray[inbe]:='-';
             SequenceTypeArray[inbe]:=False;
             SequenceArray[inbe]:=0;
             RArray[inbe]:=0;
             ValueArray[inbe]:=0;
             ValueArrayB[inbe]:=0;
             SepValueArray[inbe]:=0;
             BestFitArray[inbe]:=0;
             BestFitPhase[inbe]:=-999;
             ChoiceArray[inbe]:=0;
     end; {for inbe:=}
   PosA:=0;
   lowertrail:=0;
   PosR:=0;
   PosT:=0;
   PosB:=0;

   PhasedN:=0;
   CurrentRepresentation:=0;
   Curpair:=0;
   GCount:=0;
   ACount:=0;
   TCount:=0;
   CCount:=0;
   NCount:=0;
   XCount:=0;
   cub:=0;
   cup:=0;
   cuu:=0;
   nub:=0;
   nup:=0;
   nuu:=0;
   rb:=0;
   adb:=0;
   Dashcount:=0;
   AATTCount:=0;
   ATTACount:=0;
   SteinFormat:=False;
end; {BinClear}

Procedure SubSequenceStats;
begin
   If SubBases>0 then begin
      If SubBends>0 then CPRecord;
      SubBendFrequency:=SubBends/SubBases;
      SubPointDensity:=SubPoints/SubBases;
      Write(DisFile,SubSequenceName);
      Write(DisFile,Tab);
      Write(DisFile,FiveF:6:4);
      Write(DisFile,Tab);
      Write(DisFile,ExonF:6:4);
      Write(DisFile,Tab);
      Write(DisFile,IntronF:6:4);
      Write(DisFile,Tab);
      Write(DisFile,ThreeF:6:4);
      Write(DisFile,Tab);
      Write(DisFile,UpsF:6:4);
      Write(DisFile,Tab);
      Write(DisFile,DwnF:6:4);
      Writeln(Disfile);
   end; {if SubBases>0;}
end; {SubSequenceStats}

Procedure TotalStats;
begin
   If TotalBases>0 then begin
      SubSequenceStats;
      TotalBendFrequency:=0;
      TotalPointDensity:=0;
      If TotalBases>0 then
      TotalBendFrequency:=TotalBends/TotalBases;
      If TotalBases>0 then
      TotalPointDensity:=TotalPoints/TotalBases;
      Write(DisFile,'Whole Sequence Statistics (' + SeqFileName0 + ') Bases: ');
      Write(DisFile,TotalBases);
      Write(DisFile,' / Phased: ');
      Write(DisFile,TotalBends);
      Write(DisFile,' / Frequency: ');
      Writeln(DisFile,TotalBendFrequency);
      Write(DisFile,' / Density: ');
      Writeln(DisFile,TotalPointDensity);
      If AlliN<>0 then begin
      Writeln(DisFile,'**********************************');
Writeln(DisFile,'Whole Sequence Stats: Stein Format');
Write(DisFile,'Feature Type');
Write(DisFile,Tab);
Write(DisFile,'Phased Bases');
Write(DisFile,Tab);
Write(DisFile,'Total Bases');
Write(DisFile,Tab);
Write(DisFile,'Ratio (no char if n=0)');
Writeln(DisFile,Tab);

Write(DisFile,'Five Prime (after 300 base buffer)');
Write(DisFile,Tab);
Write(DisFile,AllfP);
Write(DisFile,Tab);
Write(DisFile,AllfN);
Write(DisFile,Tab);
If AllfN<>0 then Write(DisFile,AllfP/AllfN);
Writeln(DisFile,Tab);

Write(DisFile,'5 prime UTR (up to 200bp per gene)');
Write(DisFile,Tab);
Write(DisFile,AlluP);
Write(DisFile,Tab);
Write(DisFile,AlluN);
Write(DisFile,Tab);
If AlluN<>0 then Write(DisFile,AlluP/AlluN);
Writeln(DisFile,Tab);

Write(DisFile,'Exons');
Write(DisFile,Tab);
Write(DisFile,AlleP);
Write(DisFile,Tab);
Write(DisFile,AlleN);
Write(DisFile,Tab);
If AlleN<>0 then Write(DisFile,AlleP/AlleN);
Writeln(DisFile,Tab);

Write(DisFile,'Introns');
Write(DisFile,Tab);
Write(DisFile,AlliP);
Write(DisFile,Tab);
Write(DisFile,AlliN);
Write(DisFile,Tab);
If AlliN<>0 then Write(DisFile,AlliP/AlliN);
Writeln(DisFile,Tab);

Write(DisFile,'3 prime UTR (up to 400bp per gene)');
Write(DisFile,Tab);
Write(DisFile,AlldP);
Write(DisFile,Tab);
Write(DisFile,AlldN);
Write(DisFile,Tab);
If AlldN<>0 then Write(DisFile,AlldP/AlldN);
Writeln(DisFile,Tab);

Write(DisFile,'3 prime UTR (up to 700bp per gene)');
Write(DisFile,Tab);
Write(DisFile,AllhP);
Write(DisFile,Tab);
Write(DisFile,AllhN);
Write(DisFile,Tab);
If AllhN<>0 then Write(DisFile,AllhP/AllhN);
Writeln(DisFile,Tab);

Writeln(DisFile,'**********************************');
end; {if steinformat}
   end;
end; {TotalStats}


Procedure SubSequenceEndHandle;
var
  StillCopying: Boolean;
  inse: longint;
begin
   If PosB<PosT then PurgeArrays;
   If SubSequenceNumber>0 then SubSequenceStats;
   CPClear;
   BinClear;
   fwc1:=fwe1+fws1+100;
   SubBends:=0;
   SubBases:=0;
   SubPoints:=0;
   PosS:=0;
   PosR:=0;
   PosT:=0;
   PosB:=0;
   CurType:='f';
   OldType:='h';
   PhasedN:=0;
   BinNumber:=1;
   CommonName:='NoCommonName';
   If Curline[2]='!' then begin
      Tn1:=0;Tn2:=0;Tn3:=0;
   	  SteinFormat:=True;
      Readln(SeqFile,ExternalIndex);
      Readln(SeqFile,SubSequenceName);
      Write(SubSequenceName);
      Write(Tab);
      Readln(SeqFile,Chromosome);
      Readln(SeqFile,SubSequenceStart);
      PosW:=SubSequenceStart;
      Readln(SeqFile,Direction);
      Readln(SeqFile,CommonName);
      Write(CommonName);
      Write(Tab);
      Readln(SeqFile,SteinLength);
      Readln(SeqFile,BinType);
      FiveFiveAN:=200;FiveGeneN:=1000;GeneAGeneN:=1200;GeneThreeN:=9999999;ThreeAThreeN:=9999999;
   end
   else {i.e. if Curline[2]<>'!'}
    begin
     If Curline[2]='&' then begin
      Readln(SeqFile,ExternalIndex);
      Readln(SeqFile,SubSequenceName);
      Write(SubSequenceName);
      Write(Tab);
      Readln(SeqFile,Chromosome);
      Readln(SeqFile,SubSequenceStart);
      PosW:=SubSequenceStart;
      Readln(SeqFile,Direction);
      Readln(SeqFile,CommonName);
      Write(CommonName);
      Write(Tab);
      FiveFiveAN:=200;FiveGeneN:=1000;GeneAGeneN:=1200;GeneThreeN:=9999999;ThreeAThreeN:=9999999;
   end;
     StillCopying:=True;
     Inse:=1;
     If Length(Curline)>1 then SubSequenceName:=Curline[2];
     Inse:=3;
      While StillCopying do begin
        If (inse>length(Curline)) or (Curline[inse]=' ') then StillCopying:=False;
        If StillCopying then
        SubSequenceName:=SubSequenceName+Curline[inse];
        inse:=inse+1;
      end;
      If SubSequenceName=' ' then SubSequenceName:=SeqFileName0;
     If Length(SubSequenceName)=0 then SubSequenceName:=SeqFileName0;
     ExternalIndex:=PosW;
     BinType:='Unknown';
     Chromosome:='Unknown';
     SubSequenceStart:=PosW;
     Direction:=+1;
   end; {If Curline[2]='^'}
   SubSequenceNumber:=SubSequenceNumber+1;
   Write(Numfile,'>');
   Writeln(Numfile,SubSequenceName);
   FilterCount:=1;
end; {SubSequenceEndHandle}

Procedure DumpSource;  {This doesn't work in FPC... bummer}
var
   StillReading:Boolean;
   SourceLine:String;
begin
   Writeln(Disfile,'******************');
   Writeln(Disfile,'"' + Version + '" Source Code Follows if Available... not available in FPC: ');
   {StillReading:=True;
   SourceFileName:=Version+'.bak';
   Assign(SourceFile,SourceFileName);
   Reset(SourceFile);
   While StillReading do begin
      ReadLn(SourceFile,SourceLine);
      If EOF(SourceFile) then StillReading:=False;
      If StillReading then begin
        Write(Disfile,SourceLine);
        Writeln(DisFile);
        end; {if stillreading}
      end; {While StillReading}
   Writeln(Disfile,'******************');
   Close(SourceFile);}
end; {DumpSource}


Procedure InitHandle;
begin
   Writeln(Disfile,'Description: '+Description);
   Writeln(Disfile,'Date: '+TheDate);
   Writeln(Disfile,'Time: '+TheTime);
   Writeln(Disfile,'Sequence: '+SeqFileName0);
   Writeln(Disfile,'Sequence File Format: '+SeqFileFormat);
   Writeln(Disfile,'Sequence File Format: '+SeqFileFormat);
   Writeln(Disfile,'Analysis Program Version: ' +Version);
   Write(Disfile,'Threshold limit for phasing call: ');
   Writeln(Disfile,Threshold);
   Write(Disfile,'Random Character: ');
   Writeln(Disfile,RandomCharacter);
   Write(Disfile,'Scrambler (AA/TT <-> AT/TA switch): (Scrambled=-1, Normal=1)');
   Writeln(Disfile,Scrambler);
   Write(Disfile,'Bin Size (-1 for "feature-by-feature"): ');
   Writeln(Disfile,BinSize);
   Write(Disfile,'Number of iterations ahead the program will look for optimizing bend element: ');
   Writeln(Disfile,MoveNumber);
   Write(Disfile,'Four AA/TT Bond Value: ');
   Writeln(Disfile,FiveBondProfit);
   Write(Disfile,'Three AA/TT Bond Value: ');
   Writeln(Disfile,FourBondProfit);
   Write(Disfile,'Two AA/TT Bond Value: ');
   Writeln(Disfile,ThreeBondProfit);
   Write(Disfile,'9 Base Helical Cost: ');
   Writeln(Disfile,NineBaseLoss);
   Write(Disfile,'10 Base Helical Cost: ');
   Writeln(Disfile,TenBaseLoss);
   Write(Disfile,'11 Base Helical Cost: ');
   Writeln(Disfile,ElevenBaseLoss);
   Write(Disfile,'12 Base Helical Cost: ');
   Writeln(Disfile,TwelveBaseLoss);
   Write(Disfile,'Sequence Window Size: ');
   Writeln(Disfile,ActiveWindow);
   Writeln(Disfile);
   Writeln(Disfile,'Output Data Structure (tab delimited fields)');
   Writeln(Disfile,' 1: Common name for gene (or "NoCommonName")');
   Writeln(Disfile,' 2: SubSequence Name)');
   Writeln(Disfile,' 3: IntronF-Fraction of phased bases in intron areas');
   Writeln(Disfile,' 4: ExonF-Fraction of phased bases in exon areas');
   Writeln(Disfile,' 5: FiveF-Fraction of phased bases in 5-prime areas');
   Writeln(Disfile,' 6: ThreeF-Fraction of phased bases in 3-prime areas');
   Writeln(Disfile,' 7: IntronD-Average Phasing Score in intron areas');
   Writeln(Disfile,' 8: ExonD-Average Phasing Score in exon areas');
   Writeln(Disfile,' 9: FiveD-Average Phasing Score in 5-prime areas');
   Writeln(Disfile,' 10: ThreeD-Average Phasing Score in 3-prime areas');
   Writeln(Disfile,' 11: IntronB-Phased Bases in Introns');
   Writeln(Disfile,' 12: ExonB-Phased Bases in Exons');
   Writeln(Disfile,' 13: FiveB-Phased Bases in 5-prime');
   Writeln(Disfile,' 14: ThreeB-Phased Bases in 3-prime');
   Writeln(Disfile,' 15: IntronN-Total Bases in Introns');
   Writeln(Disfile,' 16: ExonN-Total Bases in Exons');
   Writeln(Disfile,' 17: FiveN-Total Bases in 5-prime');
   Writeln(Disfile,' 18: ThreeN-Total Bases in 3-prime');
   Writeln(Disfile,' 19: IntronP-Raw Sum of all raw Phasing Scores in Introns');
   Writeln(Disfile,' 20: ExonP');
   Writeln(Disfile,' 21: FiveP');
   Writeln(Disfile,' 22: ThreeP');
   Writeln(Disfile,' 23: IntronAF-Fraction of phased bases in Introns; First 200 nt from ATG');
   Writeln(Disfile,' 24: ExonAF');
   Writeln(Disfile,' 25: FiveAF');
   Writeln(Disfile,' 26: ThreeAF');
   Writeln(Disfile,' 27: IntronAD');
   Writeln(Disfile,' 28: ExonAD');
   Writeln(Disfile,' 29: FiveAD');
   Writeln(Disfile,' 30: ThreeAD');
   Writeln(Disfile,' 31: IntronAB');
   Writeln(Disfile,' 32: ExonAB');
   Writeln(Disfile,' 33: FiveAB');
   Writeln(Disfile,' 34: ThreeAB');
   Writeln(Disfile,' 35: IntronAN');
   Writeln(Disfile,' 36: ExonAN');
   Writeln(Disfile,' 37: FiveAN');
   Writeln(Disfile,' 38: ThreeAN');
   Writeln(Disfile,' 39: IntronAP');
   Writeln(Disfile,' 40: ExonAP');
   Writeln(Disfile,' 41: FiveAP');
   Writeln(Disfile,' 42: ThreeAP');
   Writeln(Disfile,' 43: SubSequence Number)');
   Writeln(Disfile,' 44: Bin Number)');
   Writeln(Disfile,' 45: Bin Type)');
   Writeln(Disfile,' 46: Bin Size)');
   Writeln(Disfile,' 47: Bin Start in Sub-sequence)');
   Writeln(Disfile,' 48: Bin Start in Whole sequence)');
   Writeln(Disfile,' 49: G bases)');
   Writeln(Disfile,' 50: A bases)');
   Writeln(Disfile,' 51: T bases)');
   Writeln(Disfile,' 52: C bases)');
   Writeln(Disfile,' 53: N bases)');
   Writeln(Disfile,' 54: X bases)');
   Writeln(Disfile,' 55: - bases)');
   Writeln(Disfile,' 56: AA or TT di-bases)');
   Writeln(Disfile,' 57: AT or TA di-bases)');
   Writeln(Disfile,' 58: Orientation of feature [+: same direction as sequence file, -: reversed]');
   Writeln(Disfile,' 59: Bases with phasing prediction > Threshold');
   Writeln(Disfile,' 60: Phasing Frequency');
   Writeln(Disfile,' 61: Total Value of Phasing Predictions for Bin');
   Writeln(Disfile,' 62: Phasing "Density"');
   Writeln(Disfile,' 63: cub= all capitalized bases');
   Writeln(Disfile,' 64: cup= phased capitalized bases');
   Writeln(Disfile,' 65: cuu= unphased capitailized bases');
   Writeln(Disfile,' 66: nub= all non-capitalized bases');
   Writeln(Disfile,' 67: nup= phased non-capitalized bases');
   Writeln(Disfile,' 68: nuu= unphased non-capitailized bases');
   Writeln(Disfile,' 69: adb=all determined bases (non-N) in bin');
   Writeln(Disfile,' 70: rb= all "x" bases in bin');
   Writeln(Disfile,' 71: upsF-Fraction of phased bases in ups areas');
   Writeln(Disfile,' 72: dwnF-Fraction of phased bases in dwn areas');
   Writeln(Disfile,' 73: upsD-Average Phasing Score in ups areas');
   Writeln(Disfile,' 74: dwnD-Average Phasing Score in dwn areas');
   Writeln(Disfile,' 75: upsB-Phased Bases in upss');
   Writeln(Disfile,' 76: dwnB-Phased Bases in dwns');
   Writeln(Disfile,' 77: upsN-Total Bases in upss');
   Writeln(Disfile,' 78: dwnN-Total Bases in dwns');
   Writeln(Disfile,' 79: upsP-Raw Sum of all raw Phasing Scores in upss');
   Writeln(Disfile,' 80: dwnP');
   Writeln(Disfile,' 81: External Index (position in genome)');
   {Change these definitions if the output format changes');}
   Writeln(Disfile);
   Writeln(Disfile, 'Following are the Subsequences read and evaluated by the program');
   Writeln(Disfile);
   Write(DisFile,'SequenceName');
   Write(DisFile,Tab);
   Write(DisFile,'FiveF');
   Write(DisFile,Tab);
   Write(DisFile,'ExonF');
   Write(DisFile,Tab);
  Write(DisFile,'IntronF');
  Write(DisFile,Tab);
  Write(DisFile,'ThreeF');
  Write(DisFile,Tab);
  Write(DisFile,'UpsF');
  Write(DisFile,Tab);
  Writeln(DisFile,'DwnF');
end;

Procedure UnpackRepresentation(RInt:longint);
Var idde,jdde,kdde:longint;
  begin
  jdde:=RInt;
  BaseString:='-----';
  For idde:=1 to 5 do begin
  kdde:=jdde mod 4;
  jdde:=jdde div 4;
  Case kdde of
  0:BaseString[6-idde]:='G';
  1:BaseString[6-idde]:='A';
  2:BaseString[6-idde]:='T';
  3:BaseString[6-idde]:='C';
  end; {Case}
  end; {For idde:=}
  end; {UnpackRepresentation}



Procedure Bigseq;
        var RunDone,CommentRead,KeepGoing,StillInLoop: Boolean;
        curbase:longint;
        CurChar: Char;
        Sxde:longint; {Position in file name that is being copied or changed}
        Hnde,Inde,ednI,Jnde,Knde,Lnde,Mnde,Nnde,Onde,Pnde,Qnde,Rnde,Snde,Tnde,Unde,Vnde,Wnde,Xnde,Ynde,Znde,ZndeV:longint; {Standard Indexes}
        PossibilityIndex: LongInt; {Current Possibility being evaluated}


    RepresentationIndex: LongInt; {Current representation being evaluated}
    BestScoreRecorded:Boolean;
    cvde,dvde,evde,fvde:longint;
    CurrentMove: longint; {The current Move}
begin
WriteLn('Beginning BigSeq Routine');
{Make file names for output}
{Autorun stuff: Autorun is a routine to only analyze the phased regions of worm autosomes}
AllfN:=0;AllfP:=0;AlluN:=0;AlluP:=0;AlleN:=0;AlleP:=0;AlliN:=0;AlliP:=0;AlldN:=0;AlldP:=0;AllhN:=0;AllhP:=0;
CPClear;
TPClear;
CurType:='f';
OldType:='h';
CommonName:='NoCommonName';
If SeqFileName0='Elegans' then begin
     SeqFileName:=WorkingDirectory+'x.dna';
     Autorun:=True;
     SeqFileNum:=1;
     XBase1:= 1000000;
     XBase2:=22000000;
End; {If SeqFileName=}
If SeqFileName='Autosomes' then begin
     SeqFileName:=WorkingDirectory+'i.dna';
     Autorun:=True;
     SeqFileNum:=2;
     XBase1:= 3700000;
     XBase2:=10800000;
End; {If SeqFileName=}


DataFileName:=FileNameBase+'patc';
WriteLn(DataFileName);
OFileName:=Datafilename+'.out';
DisFileName:=Datafilename+'.txt';
NumFileName:=Datafilename+'.num';
DataFileName:=Datafilename+'.dat';
SubSequenceName:=SeqFileName0;
SubSequenceNumber:=0;
WriteLn('Files about to be opened');
Assign(Seqfile,SeqFilename);
Reset(Seqfile);
Assign(Datafile,Datafilename);
Rewrite(Datafile);
Assign(Disfile,Disfilename);
Rewrite(Disfile);
Assign(Ofile,Ofilename);
Rewrite(Ofile);
Assign(NumFile,Numfilename);
Rewrite(NumFile);
Rundone:=False;
Write(Datafile,'1: Common name for gene (or "NoCommonName"');
Write(Datafile,Tab); Write(Datafile,'2: SubSequence Name');
Write(Datafile,Tab); Write(Datafile,'3: IntronF-Fraction of phased bases in intron areas');
Write(Datafile,Tab); Write(Datafile,'4: ExonF-Fraction of phased bases in exon areas');
Write(Datafile,Tab); Write(Datafile,'5: FiveF-Fraction of phased bases in 5-prime areas');
Write(Datafile,Tab); Write(Datafile,'6: ThreeF-Fraction of phased bases in 3-prime areas');
Write(Datafile,Tab); Write(Datafile,'7: IntronD-Average Phasing Score in intron areas');
Write(Datafile,Tab); Write(Datafile,'8: ExonD-Average Phasing Score in exon areas');
Write(Datafile,Tab); Write(Datafile,'9: FiveD-Average Phasing Score in 5-prime areas');
Write(Datafile,Tab); Write(Datafile,'10: ThreeD-Average Phasing Score in 3-prime areas');
Write(Datafile,Tab); Write(Datafile,'11: IntronB-Phased Bases in Introns');
Write(Datafile,Tab); Write(Datafile,'12: ExonB-Phased Bases in Exons');
Write(Datafile,Tab); Write(Datafile,'13: FiveB-Phased Bases in 5-prime');
Write(Datafile,Tab); Write(Datafile,'14: ThreeB-Phased Bases in 3-prime');
Write(Datafile,Tab); Write(Datafile,'15: IntronN-Total Bases in Introns');
Write(Datafile,Tab); Write(Datafile,'16: ExonN-Total Bases in Exons');
Write(Datafile,Tab); Write(Datafile,'17: FiveN-Total Bases in 5-prime');
Write(Datafile,Tab); Write(Datafile,'18: ThreeN-Total Bases in 3-prime');
Write(Datafile,Tab); Write(Datafile,'19: IntronP-Raw Sum of all raw Phasing Scores in Introns');
Write(Datafile,Tab); Write(Datafile,'20: ExonP');
Write(Datafile,Tab); Write(Datafile,'21: FiveP');
Write(Datafile,Tab); Write(Datafile,'22: ThreeP');
Write(Datafile,Tab); Write(Datafile,'23: IntronAF-Fraction of phased bases in Introns); First 200 nt from ATG');
Write(Datafile,Tab); Write(Datafile,'24: ExonAF');
Write(Datafile,Tab); Write(Datafile,'25: FiveAF');
Write(Datafile,Tab); Write(Datafile,'26: ThreeAF');
Write(Datafile,Tab); Write(Datafile,'27: IntronAD');
Write(Datafile,Tab); Write(Datafile,'28: ExonAD');
Write(Datafile,Tab); Write(Datafile,'29: FiveAD');
Write(Datafile,Tab); Write(Datafile,'30: ThreeAD');
Write(Datafile,Tab); Write(Datafile,'31: IntronAB');
Write(Datafile,Tab); Write(Datafile,'32: ExonAB');
Write(Datafile,Tab); Write(Datafile,'33: FiveAB');
Write(Datafile,Tab); Write(Datafile,'34: ThreeAB');
Write(Datafile,Tab); Write(Datafile,'35: IntronAN');
Write(Datafile,Tab); Write(Datafile,'36: ExonAN');
Write(Datafile,Tab); Write(Datafile,'37: FiveAN');
Write(Datafile,Tab); Write(Datafile,'38: ThreeAN');
Write(Datafile,Tab); Write(Datafile,'39: IntronAP');
Write(Datafile,Tab); Write(Datafile,'40: ExonAP');
Write(Datafile,Tab); Write(Datafile,'41: FiveAP');
Write(Datafile,Tab); Write(Datafile,'42: ThreeAP');
Write(Datafile,Tab); Write(Datafile,'43: SubSequence Number');
Write(Datafile,Tab); Write(Datafile,'44: Bin Number');
Write(Datafile,Tab); Write(Datafile,'45: Bin Type');
Write(Datafile,Tab); Write(Datafile,'46: Bin Size');
Write(Datafile,Tab); Write(Datafile,'47: Bin Start in Sub-sequence');
Write(Datafile,Tab); Write(Datafile,'48: Bin Start in Whole sequence');
Write(Datafile,Tab); Write(Datafile,'49: G bases');
Write(Datafile,Tab); Write(Datafile,'50: A bases');
Write(Datafile,Tab); Write(Datafile,'51: T bases');
Write(Datafile,Tab); Write(Datafile,'52: C bases');
Write(Datafile,Tab); Write(Datafile,'53: N bases');
Write(Datafile,Tab); Write(Datafile,'54: X bases');
Write(Datafile,Tab); Write(Datafile,'55: - bases');
Write(Datafile,Tab); Write(Datafile,'56: AA or TT di-bases');
Write(Datafile,Tab); Write(Datafile,'57: AT or TA di-bases');
Write(Datafile,Tab); Write(Datafile,'58: Orientation of feature [+: same direction as sequence file,-: reversed]');
Write(Datafile,Tab); Write(Datafile,'59: Bases with phasing prediction > Threshold');
Write(Datafile,Tab); Write(Datafile,'60: Phasing Frequency');
Write(Datafile,Tab); Write(Datafile,'61: Total Value of Phasing Predictions for Bin');
Write(Datafile,Tab); Write(Datafile,'62: Phasing "Density"');
Write(Datafile,Tab); Write(Datafile,'63: cub= all capitalized bases');
Write(Datafile,Tab); Write(Datafile,'64: cup= phased capitalized bases');
Write(Datafile,Tab); Write(Datafile,'65: cuu= unphased capitailized bases');
Write(Datafile,Tab); Write(Datafile,'66: nub= all non-capitalized bases');
Write(Datafile,Tab); Write(Datafile,'67: nup= phased non-capitalized bases');
Write(Datafile,Tab); Write(Datafile,'68: nuu= unphased non-capitailized bases');
Write(Datafile,Tab); Write(Datafile,'69: adb=all determined bases (non-N in bin');
Write(Datafile,Tab); Write(Datafile,'70: rb= all "x" bases in bin');
Write(Datafile,Tab); Write(Datafile,'71: upsF-Fraction of phased bases in ups areas');
Write(Datafile,Tab); Write(Datafile,'72: dwnF-Fraction of phased bases in dwn areas');
Write(Datafile,Tab); Write(Datafile,'73: upsD-Average Phasing Score in ups areas');
Write(Datafile,Tab); Write(Datafile,'74: dwnD-Average Phasing Score in dwn areas');
Write(Datafile,Tab); Write(Datafile,'75: upsB-Phased Bases in upss');
Write(Datafile,Tab); Write(Datafile,'76: dwnB-Phased Bases in dwns');
Write(Datafile,Tab); Write(Datafile,'77: upsN-Total Bases in upss');
Write(Datafile,Tab); Write(Datafile,'78: dwnN-Total Bases in dwns');
Write(Datafile,Tab); Write(Datafile,'79: upsP-Raw Sum of all raw Phasing Scores in upss');
Write(Datafile,Tab); Write(Datafile,'80: dwnP');
Write(Datafile,Tab); Writeln(Datafile,'81: External Index (position in genome');
PosL:=999;
PosW:=0;
PosS:=0;
PosA:=0;
lowertrail:=0;
PosB:=0;
PhasedN:=0;
PosR:=0;
PosT:=0;
TotalBends:=0;
TotalBases:=0;
TotalPoints:=0;
SubBends:=0;
SubBases:=0;
SubPoints:=0;
RetroChar:='G';
OldRetroChar:='G';
CurrentRepresentation:=0;
Curline:='>';
BinNumber:=1;
CurrentMove:=1;
For Inde:=0 to 14 do For Jnde:=0 to 4 do PhaseArray[Inde,Jnde]:=0;
For Inde:=0 to 14 do For Jnde:=0 to 16 do PhaseArray2[Inde,Jnde]:=0;
For Inde:=0 to TwoTimesActiveWindowMinusOne do begin SeparationArray[inde]:=0;  SeparationNArray[inde]:=0; SNPArray[Inde]:=0; end;
WriteLn('Setting Up Representation Index');
For RepresentationIndex:=0 to 1023 do begin
    OffBaseArray[RepresentationIndex]:=0;
    OnBaseArray[RepresentationIndex]:=0;
    Jnde:=RepresentationIndex;
    Mnde:=0;
    For Inde:=1 to 4 do begin
        Knde:=Jnde mod 4;
        Jnde:=Jnde div 4;
        Lnde:=Jnde mod 4;
        If (Scrambler=1) and ((Knde=1) or (Knde=2)) and (Lnde=Knde) then mnde:=mnde+1;
        If (Scrambler=-1) and (((Knde=1) and (Lnde=2)) or ((Knde=2) and (Lnde=1))) then mnde:=mnde+1;
    End;
    Case Mnde of
    1: begin RepresentationArray[RepresentationIndex]:=TwoBondProfit; AddArray[RepresentationIndex]:=-1; end;
    2: begin RepresentationArray[RepresentationIndex]:=ThreeBondProfit; AddArray[RepresentationIndex]:=1; end;
    3: begin RepresentationArray[RepresentationIndex]:=FourBondProfit; AddArray[RepresentationIndex]:=2; end;
    4: begin RepresentationArray[RepresentationIndex]:=FiveBondProfit; AddArray[RepresentationIndex]:=2; end;
    else
       begin RepresentationArray[RepresentationIndex]:=NoBondProfit;  AddArray[RepresentationIndex]:=-2; end;{Mnde=0}
    end {case};
end; {For representationIndex:=}

For Inde:=0 to HistogramMax do begin Histogram[Inde]:=0; LookBackArray[Inde]:=0; end;
UxA:=Histogram;FxA:=Histogram;IxA:=Histogram;ExA:=Histogram;HxA:=Histogram;DxA:=Histogram;UnA:=Histogram;FnA:=Histogram;InA:=Histogram;EnA:=Histogram;HnA:=Histogram;DnA:=Histogram;InA50:=Histogram;IxA50:=Histogram;EnA50:=Histogram;ExA50:=Histogram;
Tn1:=0;Tn2:=0;Tn3:=0;
bendsInBin:=0;

PointsInBin:=0;
   UpsB:=0;DwnB:=0;IntronB:=0;ExonB:=0;ThreeB:=0;FiveB:=0;
   UpsN:=0;DwnN:=0;IntronN:=0;ExonN:=0;ThreeN:=0;FiveN:=0;
   UpsP:=0;DwnP:=0;IntronP:=0;ExonP:=0;ThreeP:=0;FiveP:=0;
   IntronAB:=0;ExonAB:=0;ThreeAB:=0;FiveAB:=0;
   IntronAN:=0;ExonAN:=0;ThreeAN:=0;FiveAN:=0;
   IntronAP:=0;ExonAP:=0;ThreeAP:=0;FiveAP:=0;

WriteLn('Setting Up Possibility Index');

For PossibilityIndex:=1 to FourToTheMoveNumber do begin
        Jnde:=PossibilityIndex-1;
        Knde:=0;
        For Inde:=1 to Movenumber do begin
           Lnde:=Jnde mod 4;
           Jnde:=Jnde div 4;
           If inde=1 then begin
           Case Lnde of
             0: begin PossibilityArray[PossibilityIndex].MArray[Inde]:=9; PossibilityArray[PossibilityIndex].Cost[inde]:=NineBaseLoss; end;
             1: begin PossibilityArray[PossibilityIndex].MArray[Inde]:=10; PossibilityArray[PossibilityIndex].Cost[inde]:=TenBaseLoss; end;
             2: begin PossibilityArray[PossibilityIndex].MArray[Inde]:=11; PossibilityArray[PossibilityIndex].Cost[inde]:=ElevenBaseLoss; end;
             3: begin PossibilityArray[PossibilityIndex].MArray[Inde]:=12; PossibilityArray[PossibilityIndex].Cost[inde]:=TwelveBaseLoss; end;
           end; {Case}
           end; {if inde=1}
           If inde>1 then begin
           Case Lnde of
             0: begin PossibilityArray[PossibilityIndex].MArray[Inde]:=PossibilityArray[PossibilityIndex].MArray[Inde-1]+9; PossibilityArray[PossibilityIndex].Cost[inde]:=NineBaseLoss; end;
             1: begin PossibilityArray[PossibilityIndex].MArray[Inde]:=PossibilityArray[PossibilityIndex].MArray[Inde-1]+10; PossibilityArray[PossibilityIndex].Cost[inde]:=TenBaseLoss; end;
             2: begin PossibilityArray[PossibilityIndex].MArray[Inde]:=PossibilityArray[PossibilityIndex].MArray[Inde-1]+11; PossibilityArray[PossibilityIndex].Cost[inde]:=ElevenBaseLoss; end;
             3: begin PossibilityArray[PossibilityIndex].MArray[Inde]:=PossibilityArray[PossibilityIndex].MArray[Inde-1]+12; PossibilityArray[PossibilityIndex].Cost[inde]:=TwelveBaseLoss; end;
           end; {Case}
           end; {if inde>1}
        end; {For Inde}
end; {For PossibilityIndex}
ShortPossibilityArray[9]:=NineBaseLoss;
ShortPossibilityArray[10]:=TenBaseLoss;
ShortPossibilityArray[11]:=ElevenBaseLoss;
ShortPossibilityArray[12]:=TwelveBaseLoss;
ExternalIndex:=PosW;
BinType:='Unknown';
Chromosome:='Unknown';
SubSequenceStart:=PosW;
Direction:=+1;InitHandle;
BinClear;


{Reset best scores}
For inde:=1 to BestBendsRetained do begin
BestScores[inde]:=0;
BestPositions[inde]:=0;
BestStrings[inde]:='u';
end; {For inde:=1 to BestBendsRetained}
WorstBestScore:=0;
WorstBest:=1;

WriteLn('About to start loop');
           While Not(Rundone) do begin
             If PosL>length(Curline) then begin
             CommentRead:=False;
             Repeat begin
             Readln(SeqFile,Curline);
             If (length(Curline)>0) and ((Curline[1]='>') or (Curline[1]=';') or (Curline='.&') or (Curline='.!')) then begin
                If CommentRead=False then SubSequenceEndHandle; {This deals with the first comment line in a block only}
                If (Curline<>'.&') and (Curline<>'.!') then begin
                  Write(Curline);
                end; {if curline<>'.&'}
                CommentRead:=True;
              end {If (Curline[1]='>') or (Curline[1]=';') or (Curline='.&')}
              else {If (Curline[1]='>') or (Curline[1]=';') or (Curline='.&')}
               CommentRead:=False;
            end{If (Curline[1]='>') or (Curline[1]=';') or (Curline='.&')};
             If (length(Curline)=0) and EOF(SeqFile) then begin
             If Not(AutoRun) then Rundone:=True else begin
                  Close(SeqFile);
                  Case SeqFileNum of
                  1:begin SeqFileName:=WorkingDirectory+'i.dna'; SeqFileName0:='i.dna'; end;
                  2:begin SeqFileName:=WorkingDirectory+'ii.dna'; SeqFileName0:='ii.dna'; end;
                  3:begin SeqFileName:=WorkingDirectory+'iii.dna'; SeqFileName0:='iii.dna'; end;
                  4:begin SeqFileName:=WorkingDirectory+'iv.dna'; SeqFileName0:='iv.dna'; end;
                  5:begin SeqFileName:=WorkingDirectory+'v.dna'; SeqFileName0:='v.dna'; end;
                  6: Rundone:=True;
                  end; {Case}
                  SeqFileNum:=SeqFileNum+1;
                  Case SeqFileNum of
                  1:begin XBase1:= 1000000;XBase2:=22000000; end;
                  2:begin XBase1:= 3700000;XBase2:=10800000; end; {Actually the chromosome I numbers in this statement are irrelvant: set at beginning of program}
                  3:begin XBase1:= 4600000;XBase2:=11800000; end;
                  4:begin XBase1:= 4000000;XBase2:=10000000; end;
                  5:begin XBase1:= 3800000;XBase2:= 9400000; end;
                  6:begin XBase1:= 4400000;XBase2:=10400000; end;
                  end; {Case}
                  PosX:=1;
                  WriteLn('Opening File ' + SeqFileName0);
                  WriteLn('Position Start    Bases in Bin    Phased Bases    Fraction');
Assign(SeqFile,SeqFileName);
                  Reset(SeqFile);
             end {If not autorun};
             end; {If EOF(Seqfile)}
         Until (CommentRead=False) or (Rundone); {Repeat Statement}
         PosL:=1;
         end; {If PosL>length}
         Curchar:='X';    {Setting Curchar to X seems redundant but is necessary in FPC Pascal to deal with a bug if Curline is null}
         if length(Curline)>PosL-1 then Curchar:=Curline[PosL];
         PosL:=PosL+1;
         If not(Rundone) then begin
         Curbase:=5;
         Case Curchar of
         'G','g': begin Curbase:=0; end;
         'A','a': begin Curbase:=1; end;
         'T','t': begin Curbase:=2; end;
         'C','c': begin Curbase:=3; end;
         'N','n': begin Curbase:=0; end;
         'X','x': begin Curbase:=0; end;
         'R','r','Y','y','M','m','K','k','S','s','W','w': begin Curbase:=0; end;
         '/': begin
             If SteinFormat then CurType:=Curline[2];
             If SteinFormat and (Curline[2]='h') then begin GeneThreeN:=PosT; ThreeAThreeN:=PosT+200; end;
             If Not(SteinFormat) AND (PosB<posT) then begin PurgeArrays; BinClear; Curbase:=5; end;
             end;
         '-': Curbase:=5;
          {Assume 'n's are G's- this is a "neutral" base}
         end;
         If length(Curline)=0 then Curbase:=5; {This statement seems redundant but is necessary in Metrowerks Pascal to deal with a bug if Curline is null}
         If Curbase<5 then begin
             SequenceTextArray[PosA]:=CurChar;
             SequenceTypeArray[PosA]:=False;
             If (CurChar ='g') or (CurChar ='a') or (CurChar ='t') or (CurChar ='c') then lowertrail:=lowertrail+1 else lowertrail:=0;
             If lowertrail>4 then SequenceTypeArray[PosA]:=True;
             TypeArray[PosA]:=CurType;
            SequenceArray[PosA]:=Curbase;
             CurrentRepresentation:=(CurrentRepresentation*4) mod 1024 + Curbase;
             RArray[PosA]:=CurrentRepresentation;
             PosNA:=PosA-5;
             If PosNA<0 then PosNA:=PosNA+TwoTimesActiveWindow;
             ValueArrayB[PosA]:=RepresentationArray[CurrentRepresentation];
             if BalanceMode then ValueArray[PosA]:=ValueArrayB[PosA]-ValueArrayB[PosNA] else ValueArray[PosA]:=ValueArrayB[PosA];
             SepValueArray[PosA]:=AddArray[CurrentRepresentation];
             BestFitArray[PosA]:=0;
             BestFitPhase[PosA]:=-999;
             Choicearray[PosA]:=0;
             Znde:=ValueArray[PosA]; {Default Value to be bested, overall}
             ZndeV:=SepValueArray[PosA];
             If Not(AutoRun) OR NOT(CenterOff) OR ((Xbase1>PosX) or (Xbase2<PosX)) then
             {Note no begin here so centeroff just turns off the separation array in the center, not the overall phasing score}
             {Loop to handle separationArray}
             If (ThoroughSeparationArray or ((Znde>0) and (20*PhasedN>PosB) and RecordSeparationArray)) and (Not(SeparationArrayForLowerOnly) or SequenceTypeArray[PosA]) then begin {Only do separation array for bins with >5% phased residues so far}
                If (CurType='i') or Not(SteinFormat) or Not(SeparationArrayForExonJumpover) then begin
                FoundUpper:=False;
                PosNR:=5;
                StillInLoop:=True;
                If (PosNR>PosT-1) or (PosNR>TwoTimesActiveWindowMinusOne) then StillInLoop:=False;
                While StillInLoop do begin
                  RecordOn:=True;
                  If (SeparationArrayForExonJumpover or SeparationArrayNoExonJumpover) AND SteinFormat AND (TypeArray[PosNA]<>'i') then RecordOn:=False;
                  If SeparationArrayForLowerOnly and Not(SequenceTypeArray[PosNA]) then RecordOn:=False;
                  IsUpper:=False;
                  If (SequenceTextArray[PosNA]='G') or (SequenceTextArray[PosNA]='A') or (SequenceTextArray[PosNA]='T') or (SequenceTextArray[PosNA]='C') then IsUpper:=True;
                  If IsUpper then FoundUpper:=True;
                  If Not(FoundUpper) and SeparationArrayForExonJumpover then RecordOn:=False;
                  If FoundUpper and SeparationArrayNoExonJumpover then RecordOn:=False;
                  If Not(FoundUpper) AND Not(SequenceTypeArray[PosNA]) then begin RecordOn:=False; StillInLoop:=False; end;
                    If RecordOn then begin
                       SeparationArray[PosNR]:=SeparationArray[PosNR]+ZndeV*SepValueArray[PosNA];
                       SeparationNArray[PosNR]:=SeparationNArray[PosNR]+1;
                    end;
                  PosNR:=PosNR+1;
                  If PosNR<0 then PosNR:=TwoTimesActiveWindowMinusOne; {nb 7/15/02 this line seems superfluous}
                  If (PosNR>PosT-1) or (PosNR>TwoTimesActiveWindowMinusOne) then StillInLoop:=False;
                  PosNA:=PosNA-1;
                  If PosNA=-1 then PosNA:=TwoTimesActiveWindowMinusOne;
                end; {While StillInLoop}
                end; {if  (CurType='f') or Not(SteinFormat) or Not(SeparationArrayForExonJumpover)}
             end; {If ThoroughSeparationArray or...}
             Ynde:=0; {Number of turns tested}
             Mnde:=0; {How far back we are checking}
             KeepGoing:=False;
             If Znde>0 then KeepGoing:=True;
             While keepgoing do begin
               Snde:=PosA-Mnde;
               If Snde<0 then Snde:=Snde+TwoTimesActiveWindow;
               If ChoiceArray[Snde]=0 then begin
              	 Knde:=0; {Default benefit to be bested-this move}
           	    Unde:=0; {Benefit for first step}
           	    For edni:=1 to FourToTheMoveNumber do begin
           	    inde:=1+FourToTheMoveNumber-edni;
                 Lnde:=0; {Benefit for this Possibility to be modified as follows...}
                 For jnde:=1 to MoveNumber do begin
                      Xnde:=Mnde+PossibilityArray[inde].MArray[jnde];
                     Vnde:=PosA-Xnde;
                      If Vnde<0 then Vnde:=Vnde + TwoTimesActiveWindow;
                      If PosT-Xnde>-1 then Lnde:=Lnde+ValueArray[Vnde]-PossibilityArray[inde].Cost[Jnde]; {Note 12/17/04 changed > condition from 0 to -1}
           	     end {For Jnde};
                 If Lnde>Knde then begin
                   Knde:=Lnde;
                   Onde:=PossibilityArray[inde].MArray[1];
                   Wnde:=Mnde+Onde;
                   Xnde:=PosA-Wnde;
                   If Xnde<0 then Xnde:=Xnde+TwoTimesActiveWindow;
                   Unde:=-ShortPossibilityArray[PossibilityArray[inde].MArray[1]]+ValueArray[Xnde];
                 end; {If Lnde>Knde}
               end; {For Inde:=}

               If Knde=0 then begin
                 KeepGoing:=False;
                 ChoiceArray[Snde]:=-999;
               end;
               If Knde>0 then begin
                 KeepGoing:=True;
                 ChoiceArray[Snde]:=Onde;
                 Mnde:=Mnde+Onde;
                 Znde:=Znde+Unde;
               end; {If knde>0}

             end {If ChoiceArray[]=0}
               else
               begin {If ChoiceArray[]<>0}
                   If Choicearray[Snde]=-999 then Keepgoing:=False;
                   If Choicearray[Snde]>0 then begin
                       Onde:=ChoiceArray[Snde];
                       Wnde:=Mnde+Onde;
                       Xnde:=PosA-Wnde;
                       If Xnde<0 then Xnde:=Xnde+TwoTimesActiveWindow;
                       Unde:=-ShortPossibilityArray[Onde]+ValueArray[Xnde];
                       Mnde:=Mnde+Onde;
                       Znde:=Znde+Unde;
                       Keepgoing:=True;
                   end; {if coicearray[snde]>0}
               end; {else/If ChoiceArray[]<>0}
               Ynde:=Ynde+1;
               If Ynde>((Activewindow div 12)-Movenumber) then keepgoing:=false;
	               end; {while keepgoing}
               If Znde>0 then begin
                 Qnde:=PosA+1;
                 CurPhase:=0;
                 CurPeriod:=ChoiceArray[PosA];
                 If (CurPeriod<9) or (CurPeriod>12) then CurPeriod:=10;
                 For Pnde:=0 to Mnde do begin
                   Qnde:=Qnde-1;
                   If Qnde<0 then Qnde:=Qnde+TwoTimesActiveWindow;
                   If BestFitArray[Qnde]<Znde then begin
                       BestFitArray[Qnde]:=Znde;
                       BestFitPhase[Qnde]:=CurPhase;
                       end; {If BestFitArray[Qnde]<Znde}
                    CurPhase:=CurPhase+1;
                       If CurPhase=CurPeriod then begin
                          CurPhase:=0;
                          CurPeriod:=ChoiceArray[Qnde];
                          If (CurPeriod<9) or (CurPeriod>12) then Curperiod:=10;
                      end; { If CurPhase=CurPeriod}
                   end; {For pnde:=}
                   {Handle a dangling end}
                 end; { If Znde>0}
           If BestFitArray[PosA]>Threshold then begin
            Hnde:=Mnde;
           If Hnde>Histogrammax then Hnde:=HistogramMax;
           LookBackArray[Hnde]:=LookBackArray[Hnde]+1;
             cvde:=PosA-5;
            If Cvde<0 then cvde:=cvde+TwoTimesActiveWindow;
            {BestScores Stuff}
              If BestFitArray[PosA]>WorstBestScore then begin
              BestScoreRecorded:=False;
                For dvde:=1 to BestBendsRetained do begin
                If Not(BestScoreRecorded) then begin
                   If Abs(PosW-BestPositions[dvde])<50 then begin
                      If BestScores[dvde]<BestFitarray[PosA] then begin
                      BestPositions[dvde]:=PosW;
                      BestScores[dvde]:=BestFitarray[PosA];
                      BestStrings[dvde][1]:=SequenceTextArray[PosA];
                      cvde:=PosA;
                      evde:=1;
                      BestStrings[dvde]:=SequenceTextArray[cvde];
                      While Not(BestScoreRecorded) do begin
                         cvde:=cvde-1;
                      	 If cvde<0 then cvde:=TwoTimesActiveWindowMinusOne;
                         BestStrings[dvde]:=SequenceTextArray[cvde]+BestStrings[dvde];
                         evde:=evde+1;
                         If (evde>ActiveWindowMinusOne) or (evde>250) then BestScoreRecorded:=True;
                         end; {While Not(BestScoreRecorded)}
                      end; {If BestScores[dvde]<BestFitarray[PosA]}
                      BestScoreRecorded:=True;
                      end; {If Abs(PosW-BestPositions[dvde])<50}
					end; { If Not BestScoreRecorded}
					end; {For inde:=1 to BestBendsRetained}
                      If not(BestScoreRecorded) then begin
                        	BestPositions[WorstBest]:=PosW;
                        	BestScores[WorstBest]:=BestFitarray[PosA];
                            cvde:=PosA;
                            evde:=1;
                            BestStrings[WorstBest]:=SequenceTextArray[cvde];
                            While Not(BestScoreRecorded) do begin
                               cvde:=cvde-1;
                      	       If cvde<0 then cvde:=TwoTimesActiveWindowMinusOne;
                               BestStrings[WorstBest]:=SequenceTextArray[cvde]+BestStrings[WorstBest];
                               evde:=evde+1;
                               If  (evde>ActiveWindowMinusOne) or (evde>250) then BestScoreRecorded:=True;
                            end; {While Not(BestScoreRecorded)}
                      end; {If not(BestScoresRecorded}
					WorstBestScore:=99999;
					WorstBest:=1;
					For dvde:=1 to BestBendsRetained do begin
					   If BestScores[dvde]<WorstBestScore then begin
					     WorstBestScore:=BestScores[dvde];
					     WorstBest:=dvde;
					   end;
					end; {For dvde:=1 to BestBendsRetained}
              end; {If BestFitarray[PosA]>WorstBestScore}
            end; {If BestFitarray[PosA]>Threshold}

            PosA:=PosA+1;
            If PosA=TwoTimesActiveWindow then PosA:=0;
            PosT:=PosT+1;
            If PosT>ActiveWindow-1 then begin
                 PosD:=PosA-Activewindow;
                 If PosD<0 then PosD:=PosD+TwoTimesActiveWindow;
                 RecordStatus;
            end;
            PosX:=PosX+1;
         end; {If Curbase>5}
         end; {if not rundone}
   end; {While Not Rundone}
   PurgeArrays;
   TotalStats;
If SteinFormat then begin
Writeln(Ofile,'*************');
Writeln(Ofile,'Description of phasing values as a function of position in the gene, P=Phasing Fraction, N=Instances');
Write(Ofile,'Bases from beginning of sub-feature');
Write(Ofile,Tab);
Write(Ofile,'Upstream P');
Write(Ofile,Tab);
Write(Ofile,'Upstream N');
Write(Ofile,Tab);
Write(Ofile,'5prime UTR P');
Write(Ofile,Tab);
Write(Ofile,'5prime UTR N');
Write(Ofile,Tab);
Write(Ofile,'Exons P');
Write(Ofile,Tab);
Write(Ofile,'Exons N');
Write(Ofile,Tab);
Write(Ofile,'Introns P');
Write(Ofile,Tab);
Write(Ofile,'Introns N');
Write(Ofile,Tab);
Write(Ofile,'3Prime UTR P');
Write(Ofile,Tab);
Write(Ofile,'3Prime UTR N');
Write(Ofile,Tab);
Write(Ofile,'3Prime flank P');
Write(Ofile,Tab);
Write(Ofile,'3Prime flank N');
Write(Ofile,Tab);
Write(Ofile,'50 Base interval Start');
Write(Ofile,Tab);
Write(Ofile,'50 Base interval End');
Write(Ofile,Tab);
Write(Ofile,'Exons P50');
Write(Ofile,Tab);
Write(Ofile,'Exons N50');
Write(Ofile,Tab);
Write(Ofile,'Introns P50');
Write(Ofile,Tab);
Write(Ofile,'Introns N50');
Writeln(Ofile,Tab);
For Inde:=0 to HistogramMax do begin
   Write(Ofile,(10*(inde-100)));
   Write(Ofile,Tab);
   If FnA[inde]>200 then Write(Ofile,FxA[inde]/FnA[inde]);
   Write(Ofile,Tab);
   Write(Ofile,FnA[inde]);
   Write(Ofile,Tab);
   If UnA[inde]>200 then Write(Ofile,UxA[inde]/UnA[inde]);
   Write(Ofile,Tab);
   Write(Ofile,UnA[inde]);
   Write(Ofile,Tab);
   If EnA[inde]>200 then Write(Ofile,ExA[inde]/EnA[inde]);
   Write(Ofile,Tab);
   Write(Ofile,EnA[inde]);
   Write(Ofile,Tab);
   If InA[inde]>200 then Write(Ofile,IxA[inde]/InA[inde]);
   Write(Ofile,Tab);
   Write(Ofile,InA[inde]);
   Write(Ofile,Tab);
   If DnA[inde]>200 then Write(Ofile,DxA[inde]/DnA[inde]);
   Write(Ofile,Tab);
   Write(Ofile,DnA[inde]);
   Write(Ofile,Tab);
   If HnA[inde]>200 then Write(Ofile,HxA[inde]/HnA[inde]);
   Write(Ofile,Tab);
   Write(Ofile,HnA[inde]);
   Write(Ofile,Tab);
   Write(Ofile,inde*50);
   Write(Ofile,Tab);
   Write(Ofile,inde*50+49);
   Write(Ofile,Tab);
   If EnA50[inde]>1000 then Write(Ofile,ExA50[inde]/EnA50[inde]);
   Write(Ofile,Tab);
   Write(Ofile,EnA50[inde]);
   Write(Ofile,Tab);
   If InA50[inde]>1000 then Write(Ofile,IxA50[inde]/InA50[inde]);
   Write(Ofile,Tab);
   Write(Ofile,InA50[inde]);
   Writeln(Ofile,Tab);
   end; {Histogram Loop}
end; {if steinformat}

Writeln(Ofile,'*************');
Writeln(Ofile,'Bases from end of subseq   Bases in phased subseq with that location   Phased Bases with that location');
For Inde:=0 to CPMax do begin
   Write(Ofile,Inde*nTA+(nTa div 2)); {Starting point}
   Write(Ofile,Tab);
   Write(Ofile,TTA[inde]); {Number of bases with that location}
   Write(Ofile,Tab);
   Write(Ofile,TPA[inde]); {Number of bases that are bent}
   Writeln(Ofile);
   end; {Histogram Loop}

Writeln(Ofile,'*************');
Writeln(Ofile,'Value   Bases with that value   Bases with that lookback');
For Inde:=0 to HistogramMax do begin
   Write(Ofile,Inde); {Starting point}
   Write(Ofile,Tab);
   Write(Ofile,Histogram[inde]); {Number of bases with that bend value}
   Write(Ofile,Tab);
   Write(Ofile,LookbackArray[inde]); {Number of bases with that bend value}
   Write(Ofile,Tab);
   Writeln(Ofile);
   end; {Histogram Loop}

If RecordSeparationArray or ThoroughSeparationArray then begin
Writeln(Ofile,'*************');
Writeln(Ofile,'Separation Array');
Write(Ofile,'SeparationValue(bases)');
Write(Ofile,Tab);
Write(Ofile,'Coincidences value for that separation');
Write(Ofile,Tab);
Write(Ofile,'Opportunities');
Write(Ofile,Tab);
Write(Ofile,'Coincidence value per opportunity');
Writeln(Ofile,Tab);
For Inde:=0 to TwoTimesActiveWindowMinusOne do begin
   Write(Ofile,Inde); {Starting point}
   Write(Ofile,Tab);
   Write(Ofile,SeparationArray[inde]); {Coincidence value for that separation value}
   Write(Ofile,Tab);
   Write(Ofile,SeparationNArray[inde]); {Number of bases tested}
   Write(Ofile,Tab);
   If SeparationNArray[inde]>0 then Write(Ofile,SeparationArray[inde]/SeparationNArray[inde]); {Number of bases with that separation value}
   Write(Ofile,Tab);
   Writeln(Ofile);
   end; {Histogram Loop}
   end; {If RecordSeparationArray}

If RecordSNPArray then begin
Writeln(Ofile,'*************');
Writeln(Ofile,'Position in SNP data element (bases)   Total value of indices with that position');
For Inde:=0 to TwoTimesActiveWindowMinusOne do begin
   Write(Ofile,Inde); {Starting point}
   Write(Ofile,Tab);
   Write(Ofile,SNPArray[inde]); {Number of bases with that separation value}
   Write(Ofile,Tab);
   Writeln(Ofile);
   end; {Histogram Loop}
   end; {If RecordSeparationArray}

Writeln(Ofile,'*************');
Writeln(Ofile,'Base composition in phased (lines 1-14) and unphased (line15) regions:  Phasing ...  G ... A ... T ... C ... N');
For Inde:=0 to 14 do begin
   Write(Ofile,Inde); {Phasing}
   Write(Ofile,Tab);
   For jnde:=0 to 4 do begin
      Write(Ofile,PhaseArray[inde,jnde]); {Number of bases with that phasing value}
      Write(Ofile,Tab);
   end; {jnde Loop}
   Write(Ofile,Tab);
   Writeln(Ofile);
   end; {inde Loop}

Writeln(Ofile,'*************');
Writeln(Ofile,'Dinucleotide composition in phased (lines 1-14) and unphased (line15) regions:  Phasing . GG . GA . GT . GC . AG . AA . AT . AC . TG . TA . TT . TC . CG . CA . CT . CC . Total');
For Inde:=0 to 14 do begin
   Write(Ofile,Inde); {Phasing}
   Write(Ofile,Tab);
   For jnde:=0 to 16 do begin
      Write(Ofile,PhaseArray2[inde,jnde]); {Number of bases with that phasing value}
      Write(Ofile,Tab);
   end; {jnde Loop}
   Write(Ofile,Tab);
   Writeln(Ofile);
   end; {inde Loop}

Writeln(Ofile,'*************');
Writeln(Ofile,'Five Base Sequence   Frequency On Face of Helix   Frequency on rear of Helix');
For Inde:=0 to 1023 do begin

   Write(Ofile,Inde); {Representation Value}
   Write(Ofile,Tab);
   UnpackRepresentation(Inde);
   Write(Ofile,BaseString); {Sequence}
   Write(Ofile,Tab);
   Write(Ofile,OnBaseArray[inde]); {Number of interstitial sequences with that value}
   Write(Ofile,Tab);
   Write(Ofile,OffBaseArray[inde]); {Number of interstitial sequences with that value}
   Writeln(Ofile);
   end; {representation Loop}

Writeln(Ofile);
Writeln(Ofile,'*************');
Writeln(Ofile,'Most closely matching sequences in dataset follow- end of sequence is closest match');

For inde:=1 to BestBendsRetained do begin
    If BestPositions[inde]<>0 then begin
     Write(Ofile,'Approximate numerical position of sequence=');
     Write(Ofile,BestPositions[inde]);
     Writeln(Ofile);
     Write(Ofile,'Score of position=');
     Write(Ofile,BestScores[inde]);
     WriteLn(Ofile);
     Lnde:=1;
     Mnde:=1;
     While Lnde<ActiveWindow+1 do begin
     If length(BestStrings[inde])>Lnde-1 then Write(Ofile,BestStrings[inde][Lnde]);
     If Mnde=10 then begin
         Writeln(Ofile);
         Mnde:=0;
         end; {if mnde=10}
     Mnde:=Mnde+1;
     Lnde:=Lnde+1;
     end; {While Lnde<Activewindow-1}
     Writeln(Ofile);
	 Writeln(Ofile);
	 end; {If BestPositions[inde]<>0}
	 end; {For inde:=1 to BestBendsRetained}
Writeln(Ofile,'*************');
DumpSource;
Close(DisFile);
Close(SeqFile);
Close(DataFile);
Close(NumFile);
Close(OFile);
end; {BigSeq}

BEGIN {Main Program}
{Randomize;}
WorkingDirectory:=ParamStr(0);
While (length(WorkingDirectory)>0) and (WorkingDirectory[length(WorkingDirectory)]<>'\') and (WorkingDirectory[length(WorkingDirectory)]<>'/') do Delete(WorkingDirectory,length(WorkingDirectory),1);

Tab:= Chr(9);
ln2:=0.69314718056;
InputFileName:=WorkingDirectory+InputFileBase;
Assign(InputFile,InputFileName);
Reset(InputFile);
Readln(InputFile,SeqFileName);SeqFileName0:=SeqFileName;
if (Pos('/',SeqFileName)=0) and (Pos(':',SeqFileName)=0) then SeqFileName:=WorkingDirectory+SeqFileName;

UserStillGoing:=True;
If EOF(InputFile) then begin
   Write('Sequence File Name [e.g. ***.seq {"Autosomes": C. elegans Autosomes, "Elegans": C. elegans genome} ]? : (q=quit)');
   Readln(SeqFileName); SeqFileName0:=SeqFileName;
     If SeqFileName='q' then UserStillGoing:=False;
   if (Pos('/',SeqFileName)=0) and (Pos(':',SeqFileName)=0) then SeqFileName:=WorkingDirectory+SeqFileName;
     While UserStillGoing do begin
      Write('Output File Prefix [e.g. "OutFile"]? : ');
      Readln(FileNameBase);
     if (Pos('/',FileNameBase)=0) and (Pos(':',FileNameBase)=0) then FileNameBase:=WorkingDirectory+FileNameBase;

      Write('Bin Size?');
      Readln(BinSize);
      Write('Score Threshold for phasing (e.g., 100)? : ');
      Readln(Threshold);
      Write('Random Number Seed? ("0" for no randomization) : ');
      Readln(RandomCharacter);
      Write('Scramble the sequence ("-1": Sequence will be scrambled, "+1": Normal : ');
      Readln(Scrambler);
      Write('Describe this experiment: ');
      Readln(Description);
      Write('Todays Date? : ');
      Readln(TheDate);
      Write('What Time is it? : ');
      Readln(TheTime);
      Write('Filter window start n bases upstream of phasing (e.g. "150")');
      Readln(fws1);
      Write('Filter window end n bases downstream of phasing (e.g. "150")');
      Readln(fwe1);
      Writeln('Other Options?');
      Writeln;
      Writeln('" C " : Center Off [C. elegans]');
      Writeln('" S " : Separation array (limit to phased regions)');
      Writeln('" T " : Separation array (Thorough-all regions)');
      Writeln('" R " : Reflect CP Array');
      Writeln('" N " : Numerical Array');
      Writeln('" H " : SuperSensitive');
      Writeln('" X " : SNP Array');
      Writeln('" L " : Only include lower case for separation array');
      Writeln('" J " : Only include lower case sequences separated by at least one upper case for separation array');
      Writeln('" K " : Exclude lower case sequences separated by at least one upper case for separation array');
      Writeln('" F " : Filter Array: Generate .Num File filtered to contain ns exept where sequence has above-threshold phasing');
      Writeln('" B " : Balance Mode: Balance possible high score in AA/TT-rich regions by subtracting off scores on opposite side of helix');
      Writeln;
      {Handle ParameterList}
      Readln(ParamaterList);
      RecordSeparationArray:=False;
      RecordSNPArray:=False;
      ThoroughSeparationArray:=False;
      SeparationArrayForLowerOnly:=False;
      SeparationArrayForExonJumpover:=False;
      SeparationArrayNoExonJumpover:=False;
      CenterOff:=False;
      ReflectCPArray:=False;
      NumericalArray:=False;
      SuperSensitive:=False;
      FilterArray:=False;
      BalanceMode:=False;
      For Iade:=1 to length(ParamaterList) do begin
         Case ParamaterList[Iade] of
         's','S': RecordSeparationArray:=True;
         'x','X': RecordSNPArray:=True;
         'c','C': CenterOff:=True;
         'r','R': ReflectCPArray:=True;
         'n','N': NumericalArray:=True;
         'h','H': SuperSensitive:=True;
         't','T': ThoroughSeparationArray:=True;
         'l','L': SeparationArrayForLowerOnly:=True;
         'j','J': SeparationArrayForExonJumpover:=True;
         'k','K': SeparationArrayNoExonJumpover:=True;
         'f','F': FilterArray:=True;
         'b','B': BalanceMode:=True;

         end; {case}
         end; {For Iade:=}
      Writeln('Ready to Go');
      Bigseq;
      Writeln('Finished ' + SeqFileName0);
      Writeln('***********');
      Write('Sequence File Name [e.g. ***.seq]? : (q=quit)');
      Readln(SeqFileName);SeqFileName0:=SeqFileName;
      If SeqFileName0='q' then UserStillGoing:=False;
      if (Pos('/',SeqFileName)=0) and (Pos(':',SeqFileName)=0) then SeqFileName:=WorkingDirectory+SeqFileName;

   end; {While UserStillGoing}
Close(InputFile);
end else begin {if EOF}
    If (SeqFileName0='q') then UserStillGoing:=False;
     While UserStillGoing do begin
      Write('Output File Prefix [e.g. "OutFile"]? : ');
      Readln(InputFile,FileNameBase);
      if (Pos('/',FileNameBase)=0) and (Pos(':',FileNameBase)=0) then FileNameBase:=WorkingDirectory+FileNameBase;

      Writeln(FileNameBase);
      Write('Bin Size?');
      Readln(InputFile,BinSize);
      Writeln(BinSize);
      Write('Score Threshold for phasing (e.g., 100)? : ');
      Readln(InputFile,Threshold);
      Writeln(Threshold);
      Write('Random Number Seed? ("0" for no randomization) : ');
      Readln(InputFile,RandomCharacter);
      Writeln(RandomCharacter);
      Write('Scramble the sequence ("-1": Sequence will be scrambled, "+1": Normal : ');
      Readln(InputFile,Scrambler);
      Writeln(Scrambler);
      Write('Describe this experiment: ');
      Readln(InputFile,Description);
      Writeln(Description);
      Write('Todays Date? : ');
      Readln(InputFile,TheDate);
      Writeln(TheDate);
      Write('What Time is it? : ');
      Readln(InputFile,TheTime);
      Writeln(TheTime);
      Write('Filter window start n bases upstream of phasing (e.g. "+150")');
      Readln(InputFile,fws1);
      Writeln(fws1);
      Write('Filter window start n bases downstream of phasing (e.g. "+150")');
      Readln(InputFile,fwe1);
      Writeln(fwe1);
      {Handle ParameterList}
      Writeln('Other Options?');
      Writeln;
      Writeln('" C " : Center Off [C. elegans]');
      Writeln('" S " : Separation array (limit to phased regions)');
      Writeln('" T " : Separation array (Thorough-all regions)');
      Writeln('" R " : Reflect CP Array');
      Writeln('" N " : Numerical Array');
      Writeln('" H " : SuperSensitive');
      Writeln('" X " : SNP Array');
      Writeln('" L " : Only include lower case for separation array');
      Writeln('" J " : Only include lower case sequences separated by at least one upper case for separation array');
      Writeln('" K " : Exclude lower case sequences separated by at least one upper case for separation array');
      Writeln('" F " : Filter Array: Generate .Num File filtered to contain ns exept where sequence has above-threshold phasing');
      Writeln('" B " : Balance Mode: Balance possible high score in AA/TT-rich regions by subtracting off scores on opposite side of helix');
      Writeln;
      Readln(InputFile,ParamaterList);
      RecordSeparationArray:=False;
      RecordSNPArray:=False;
      ThoroughSeparationArray:=False;
      ReflectCPArray:=False;
      FilterArray:=False;
      CenterOff:=False;
      For Iade:=1 to length(ParamaterList) do begin
         Case ParamaterList[Iade] of
         's','S': RecordSeparationArray:=True;
         'x','X': RecordSNPArray:=True;
         't','T': ThoroughSeparationArray:=True;
         'c','C': CenterOff:=True;
         'r','R': ReflectCPArray:=True;
         'n','N': NumericalArray:=True;
         'h','H': SuperSensitive:=True;
         'l','L': SeparationArrayForLowerOnly:=True;
         'j','J': SeparationArrayForExonJumpover:=True;
         'k','K': SeparationArrayNoExonJumpover:=True;
         'f','F': FilterArray:=True;
         'b','B': BalanceMode:=True;
        end; {case}
         end; {For Iade:=}
      Writeln('******************');
      Writeln('Ready to Go');
      Bigseq;
      Writeln('Finished ' + SeqFileName);
      Writeln('*******************');
      Write('Sequence File Name [e.g. ***.seq]? : (q=quit)');
      Readln(InputFile,SeqFileName); SeqFileName0:=SeqFileName;
      If (SeqFileName0='q') Or EOF(InputFile) then UserStillGoing:=False;
      if (Pos('/',SeqFileName)=0) and (Pos(':',SeqFileName)=0) then SeqFileName:=WorkingDirectory+SeqFileName;
      end; {While UserStillGoing}
end;
WriteLn(' ');
WriteLn('All Done');
End.

{Typical Structure of Batch File 'PATC.bat'
Line 1: SequenceFileName
Line 2: OutputFileBase
Line 3: Bin Size
Line 4: Threshold for phasing
Line 5: Scramble (+1 for not, -1 for scrambled)
Line 6: Description of experiment (256 char max)
Line 7: Date
Line 8: Time
Line 9: FilterWindowStart
Line 10: FilterWindowEnd
Line 11: Options: "C":Center Off [C. elegans], "S", Separation array, "R":Reflect CP Array
Line 12: q (quit) or start new directive with SequenceFileName, etc

Can be catenatted, with the last line a "q" to indicate quit.}
