function shiSpmMultiCond_PresdLog(Logfile,TrInfo,Trigger,CondInfo,OutPath,Outfile,existSubjectID)

% reads log files produced by NBS Presentation and converts them to multiple-condition .mat files for SPM
% 
% shiSpmMultiCond_PresdLog(Logfile,TrInfo,Trigger,CondInfo,OutPath,Outfile)
% shiSpmMultiCond_PresdLog(Logfile,TrInfo,Trigger,CondInfo,OutPath,Outfile,existSubjectID)
% 
%   Logfile           - string or cell array of strings of input .log file
%                       names
%   TrInfo .Tr        - double, length of TR in second
%          .TrDiscd   - int, number of TRs discarded at the beginning
%   Trigger           - string of trigger code
%   CondInfo.CondCode - N*1 cell array of strings of codes (which can be
%                       regular expressions), N being the number of
%                       conditions
%           .EvtType  - 1*K cell array of strings of acceptable event types
%                       e.g. {'Picture','Video'}
%           .Dur      - 1*N vector of doubles, duration of each condition
%                       in TR
%           .CondName - 1*N cell array of condition names
%   OutPath           - string, path of output
%   Outfile           - string or cell array of strings of output .mat file
%                       names
%   existSubjectID    - assign 1 if the first column in the logfiles is
%                       subject ID (default = 0)
% 
%    ###########
% by Zhenhao Shi @ 2014-12-30
%    ###########
% 

if nargin < 7
    existSubjectID = 0;
end

if ischar(Logfile)
    Logfile = {Logfile};
end
if ischar(Outfile)
    Outfile = {Outfile};
end
if ischar(CondInfo.EvtType)
    CondInfo.EvtType = {CondInfo.EvtType};
end

if length(Logfile) ~= length(Outfile)
    error('number of input .log files and output .mat files not matched');
end

Tr = TrInfo.Tr;
n_DiscardedTr = TrInfo.TrDiscd;
ConditionCode = CondInfo.CondCode;
EventType = CondInfo.EvtType;
ConditionDuration = CondInfo.Dur;
Condition = CondInfo.CondName;

cwd = pwd;

if ~isdir(OutPath)
    mkdir(OutPath);
end
OutPath = shiFullFileName(OutPath);
OutPath = OutPath{1};

cd(OutPath);

n_Code = zeros(length(Condition),1);
for j = 1:length(Condition)
    n_Code(j) = length(ConditionCode{j});
end;

Vector = cell(length(Logfile),1);




    for run = 1:length(Logfile)

        fid = fopen(Logfile{run},'rt');
        if existSubjectID
            Log = textscan(fid,'%*s%*s%s%s%f%*[^\n]','delimiter','\t','headerlines',5);
        else
            Log = textscan(fid,'%*s%s%s%f%*[^\n]','delimiter','\t','headerlines',5);
        end
        LogType = Log{1};
        LogCode = Log{2};
        LogTime = Log{3};

        for m = 1:length(LogCode)
            if strcmp(LogCode{m},Trigger) && strcmp(LogType{m},'Response');
                break;
            end
        end;
        StartTime = LogTime(m);

        isStimulus = shiStrIncl(LogType,EventType);
%         isStimulus = find(isStimulus);
%         for typ = 1:length(EventType)
%             isStimulusTemp = strcmp(EventType{typ},LogType);
%             for m = 1:length(isStimulusTemp)
%                 if isStimulusTemp(m)
%                     isStimulus = [isStimulus,m];
%                 else
%                     continue;
%                 end
%             end;
%         end;
%         isStimulus = sort(isStimulus);
        StimulusCode = LogCode(isStimulus);
        StimulusTime = LogTime(isStimulus);

        for con = 1:length(Condition)

            isCondition = shiStrIncl(StimulusCode,ConditionCode{con});
%             isCondition = find(isCondition);
%             for p = 1:n_Code(con)
%                 isConditionTemp = strcmp(ConditionCode{con}{p},StimulusCode);
%                 for m = 1:length(isConditionTemp)
%                     if isConditionTemp(m)
%                         isCondition = [isCondition,m];
%                     else
%                         continue;
%                     end
%                 end;
%             end;
%             isCondition = sort(isCondition);

            Vector{run}{con}=(StimulusTime(isCondition)'-StartTime)/(Tr*10000)-n_DiscardedTr;

         end;

        fclose(fid);
        clear('Log','LogType','LogCode','LogTime','StartTime','isStimulus','isCondition','StimulusCode','StimulusTime');

    end;



fid = fopen(sprintf('Information_%s.txt',shiTime),'at');

for i = 1:length(Condition)
    fprintf(fid,'%s%s%s%d\n',['Condition',num2str(i,'%.2d'),': '],Condition{i},'; Duration=',ConditionDuration(i));
end;

fprintf(fid,'\n%s','Number of trials for each condition per log:');

    for run = 1:length(Logfile)
        fprintf(fid,'\n%s\n',Logfile{run});
        for con = 1:length(Condition)
            fprintf(fid,' %d',length(Vector{run}{con}));
        end;
    end;

fprintf(fid,'\n\n\n');
fclose(fid);


    for run = 1:length(Logfile)
        names = cell(1,length(Condition));
        onsets = cell(1,length(Condition));
        durations = cell(1,length(Condition));
        for con = 1:length(Condition)
            names{con} = Condition{con};
            onsets{con} = Vector{run}{con};
            durations{con} = ConditionDuration(con);
        end;
        save(Outfile{run},'names','onsets','durations');
        clear('names','onsets','durations');
    end;


cd(cwd);
