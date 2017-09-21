function [str, out] = sw_example(str,path,figname,recalc)
% executes example code given in a cell of strings
%
% Lines starting with >> will be executed
% Lines starting with >>> will be executed, but not shown
% Line containing snapnow causes to save the image of the last figure

keep = true(1,numel(str));

out = {};

warn0 = warning;
warning off
cc = onCleanup(@() warning(warn0));

% find code blocks
idx = 1;

for ii = 1:numel(str)
    temp = strtrim(str{ii});
    if ~isempty(regexp(temp,'snapnow','once'))
        % save figure
        fignamei = [figname '_' num2str(idx) '.png'];
        
        if recalc    
            idx = idx+1;
            set(gcf,'Color',[241 240 240]/255)
            drawnow;
            ihc = get(gcf,'InvertHardCopy');
            set(gcf,'InvertHardCopy','off');
            print([path filesep fignamei],'-dpng','-r144');
            set(gcf,'InvertHardCopy',ihc);
        end
        caption = str{ii-1};
        caption = strtrim(caption(caption~='>'));
        str{ii} = ['```' newline ' ' newline '{% include image.html file="' fignamei '" alt="' caption '" max-width="500" %}' newline newline '```matlab' newline];
        keep(ii) = true;
    elseif numel(temp)>3 && strcmp(temp(1:3),'>>>')
        if recalc
            evalc(str{ii}(4:end));
        end
        keep(ii) = false;
    elseif numel(temp)>2 && strcmp(temp(1:2),'>>')
        if recalc
            out0 = evalc(str{ii}(3:end));
        else
            out0 = '';
        end
        if ~isempty(out0)
            out{end+1} = out0; %#ok<AGROW>
        end
        str{ii} = str{ii}(find(str{ii}=='>',1)+2:end);
    end
end

% close all opened figures
close('all')
str = str(keep);
str = strsplit(sprintf('%s\n',str{:}),newline)';
str = str(1:end-1);


end