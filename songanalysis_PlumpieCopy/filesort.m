function [FileTypes] = filesort(FileName, FileExt, StartNum, EndNum, string)

fileNo = EndNum - StartNum;
BOSList = 'BOS: ';
rBOSList = 'rBOS: ';
              
BOS_No = 0;
rBOS_No = 0;

for i = 0:fileNo,
    if (StartNum+i < 10);
       RecFileName = strcat(FileName,'.00',num2str(StartNum+i),FileExt);
    else
       if (StartNum+i < 100);
       RecFileName = strcat(FileName,'.0',num2str(StartNum+i),FileExt);
        else
       RecFileName = strcat(FileName,'.',num2str(StartNum+i),FileExt);
        end
    end

    if ~(exist(RecFileName,'file'))
        continue;
    end
    
    fid = fopen(RecFileName,'r');
%   disp(RecFileName); %
    Number = 0;
    index = 0;
    while (index < 5)
       tline = fgetl(fid);
        if (feof(fid)),
            break;
        end
        
        if ((strfind(tline,string) > 0))
 %             disp('forward')%
              Number = (StartNum + i);
              BOS_No = BOS_No + 1;
              index = 20;
        end
%         if (strfind(tline,'greenorange-07232007.386.cbin.rsong') > 0)
% %              disp('reverse')%
%               Number = (StartNum + i);
%               rBOS_No = rBOS_No + 1;
%               index = 10;
%         end

%disp(index)%

        if (index > 10)
%             disp(Number)%
              if (Number < 10);
                  NumString = strcat('00',num2str(Number));
              else
                 if (Number < 100);
                  NumString = strcat('0',num2str(Number));          
                 else
                   NumString = num2str(Number);
                 end
%       disp(NumString)%
              end
            FileTypes.BOS(BOS_No,:) = Number;
            BOSList =strcat(BOSList, NumString,',');
        end
          
          if (index < 19),
             if (index > 9)
%             disp(Number)%
              if (Number < 10);
                  NumString = strcat('00',num2str(Number));
              else
                 if (StartNum+i < 100);
                  NumString = strcat('0',num2str(Number));          
                 else
                   NumString = num2str(Number);
                 end
%        disp(NumString)%
              end
%             FileTypes.RevBOS(rBOS_No,:) = Number;
%             rBOSList =strcat(rBOSList, NumString,',');              
             end
          end
    end
end
disp (BOSList)
% disp (rBOSList)
end


   
%    Stimulus: greenorange-07232007.386.song %
    
% BOS = %
% rBOS = %