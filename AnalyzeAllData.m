
% Run all sessions and create 'AnalyzedData.mat' containing all data

function AnalyzeAllData(Filename)
    
    RunAllSessionsTimeVsDistance
    load('Summary.mat');
    clear R;
    Cnt = 0;
    Mx = size(R,2)+2;
    for i=1:size(R1,1)
        if R1(i,1)~=0
            Cnt=Cnt+1;               
            R(Cnt,1)=1;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R1(i,j-2);
            end
        end
    end
    for i=1:size(R2,1)
        if R2(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=2;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R2(i,j-2);
            end
        end
    end
    for i=1:size(R3,1)
        if R3(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=3;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R3(i,j-2);
            end
        end
    end
    for i=1:size(R4,1)
        if R4(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=4;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R4(i,j-2);
            end
        end
    end
    for i=1:size(R5,1)
        if R5(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=5;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R5(i,j-2);
            end
        end
    end
    for i=1:size(R6,1)
        if R6(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=6;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R6(i,j-2);
            end
        end
    end
    for i=1:size(R7,1)
        if R7(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=7;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R7(i,j-2);
            end
        end
    end
    for i=1:size(R8,1)
        if R8(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=8;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R8(i,j-2);
            end
        end
    end
    for i=1:size(R9,1)
        if R9(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=9;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R9(i,j-2);
            end
        end
    end
    for i=1:size(R10,1)
        if R10(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=10;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R10(i,j-2);
            end
        end
    end
    for i=1:size(R11,1)
        if R11(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=11;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R11(i,j-2);
            end
        end
    end
    for i=1:size(R12,1)
        if R12(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=12;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R12(i,j-2);
            end
        end
    end
    for i=1:size(R13,1)
        if R13(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=13;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R13(i,j-2);
            end
        end
    end
    for i=1:size(R14,1)
        if R14(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=14;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R14(i,j-2);
            end
        end
    end
    for i=1:size(R15,1)
        if R15(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=15;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R15(i,j-2);
            end
        end
    end
    for i=1:size(R16,1)
        if R16(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=16;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R16(i,j-2);
            end
        end
    end
    for i=1:size(R17,1)
        if R17(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=17;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R17(i,j-2);
            end      
        end
    end
    for i=1:size(R18,1)
        if R18(i,1)~=0
            Cnt=Cnt+1;
            R(Cnt,1)=18;
            R(Cnt,2)=i;
            for j=3:Mx
                R(Cnt,j)=R18(i,j-2);
            end
        end
    end
    
    save('AnalyzedData.mat','R');
    AnalyzePopulation('AnalyzedData');
end