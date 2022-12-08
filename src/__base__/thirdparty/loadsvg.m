%% load svg

% svgname: filename you want to convert
% fine: fineness of svg path. The smaller the number, the higher the detail.
% 0<fine<1
% graph: input true or false. if graph is true display plot graph.

function R=loadsvg(svgname,fine,graph)

svgtext=fileread(svgname);

%% Prepare to trasnsform svg string to array
P_status='[MmCcLlSsVvHhZz]';
P_path_str=extractBetween(svgtext,'<path','>');
P_path_str=regexprep(P_path_str,'id="','ID="');
P_d_str=extractBetween(P_path_str,'d="','"');
P_d_str=regexprep(P_d_str,',',' ');
P_d_str=regexprep(P_d_str,'-',' -');
% % P_line_str=regexprep(P_line_str,'(?<!-)[.]',' .');
P_d_str=regexprep(P_d_str,P_status,' $0 ');
% P_line_str=strip(P_line_str);
P_d_str=regexprep(P_d_str,'\s*',' ');
% P_line_str=regexprep(P_line_str,',,',',');

num_P_cell=length(P_d_str);

P_line=[];

% .123.456.678みたいな数字を分ける。
% 配列のセル配列に変換する。
for id_cell=1:num_P_cell
    
    P_str_load=P_d_str{id_cell};
    num_dot=0;
    save_sp=false;
    num_str=1;
    
    while 1
        str_load=P_str_load(num_str);
        
        if str_load==' '
            save_sp=true;
            num_dot=0;
        end
        
        if str_load=='.'
            num_dot=num_dot+1;
            
            if save_sp
                save_sp=false;
            end
            
            if save_sp==0 && num_dot>1
                P_str_load=[P_str_load(1:num_str-1),' .',P_str_load(num_str+1:end)];
                num_str=num_str+1;
            end
        end
        
        if num_str>=strlength(P_str_load)
            P_line_dum=split(P_str_load,' ');
            Emp_cell=find(cellfun('isempty',P_line_dum));
            P_line_dum(Emp_cell)=[];
            P_line{id_cell}=P_line_dum;
            break
        end
        
        num_str=num_str+1;
    end
end

%% Convert svg text to coordinates
% R=
% col.1 | col.2 | col.3...
% m     |m      |m     ...
% 2.13  |423    |78    ...

% fine=0.2; %Set fineness of svg line. Smaller value means more fine. (0<fine<1)

for id_line=1:num_P_cell
    P_line_inuse=P_line{id_line};
    
    % Set initial point
    p_Origin=[0,0];
    p_por_origin=p_Origin;
    
    R{id_line}=[];
    
    % Set initial reference point
    pref=p_Origin;
    
    % Set initial order. I don't know svg not to begin 'M' or 'm'.
    sta=[];
    
    % Set initial num_row
%     num_row=4;
    num_row=1;
    
    while num_row<=length(P_line_inuse)
        switch P_line_inuse{num_row}
        case {'M','m','C','c','L','l','S','s','V','v','H','h'}
            sta(end+1)=P_line_inuse{num_row};
            num_row=num_row+1;
        otherwise
            if sta(end)=='M'
                if P_line_inuse{num_row-1}=='M'
                    p0=[str2double(P_line_inuse{num_row}),str2double(P_line_inuse{num_row+1})];
                    
                    p_por_origin=p0; 
                    pref=p0;
                else        
                    p0=pref;

                    p1=[str2double(P_line_inuse{num_row}),str2double(P_line_inuse{num_row+1})];

                    l=svgl(fine,p0,p1);

                    pref=p1;

                    R{id_line}=[R{id_line};l];
                end
                
                num_row=num_row+2;
                
            elseif P_line_inuse{num_row}=='Z'
                p0=pref;

                p1=p_por_origin;

                l=svgl(fine,p0,p1);

                pref=p1;

                R{id_line}=[R{id_line};l];
                
                num_row=num_row+1;
                
                sta(end+1)='Z';
                
            elseif P_line_inuse{num_row}=='z'
                p0=pref;

                p1=p_por_origin;

                l=svgl(fine,p0,p1);

                pref=p1;

                R{id_line}=[R{id_line};l];
                
                num_row=num_row+1;
                
                sta(end+1)='z';
                
            elseif sta(end)=='m'
                if P_line_inuse{num_row-1}=='m'
                    p0=[str2double(P_line_inuse{num_row}),str2double(P_line_inuse{num_row+1})]+pref;
                    p_por_origin=p0;
                    pref=p0;
                else
                    p0=pref;

                    p1=[str2double(P_line_inuse{num_row}),str2double(P_line_inuse{num_row+1})]+pref;

                    l=svgl(fine,p0,p1);

                    pref=p1;

                    R{id_line}=[R{id_line};l];
                    
%                     p_por_origin=p1;
                end
                
                num_row=num_row+2;
                
            elseif sta(end)=='L'
                p0=pref;

                p1=[str2double(P_line_inuse{num_row}),str2double(P_line_inuse{num_row+1})];

                l=svgl(fine,p0,p1);

                pref=p1;

                R{id_line}=[R{id_line};l];
                
                num_row=num_row+2;
                                
            elseif sta(end)=='l'
                p0=pref;

                p1=[str2double(P_line_inuse{num_row}),str2double(P_line_inuse{num_row+1})]+pref;

                l=svgl(fine,p0,p1);

                pref=p1;

                R{id_line}=[R{id_line};l];
                
                num_row=num_row+2;
                
            elseif sta(end)=='V'
                p0=pref;

                p1=[pref(1),str2double(P_line_inuse{num_row})];

                l=svgl(fine,p0,p1);

                pref=p1;

                R{id_line}=[R{id_line};l];
                
                num_row=num_row+1;
                
            elseif sta(end)=='v'
                p0=pref;

                p1=[0,str2double(P_line_inuse{num_row})]+pref;

                l=svgl(fine,p0,p1);

                pref=p1;

                R{id_line}=[R{id_line};l];
                
                num_row=num_row+1;
                
            elseif sta(end)=='H'
                p0=pref;

                p1=[str2double(P_line_inuse{num_row}),pref(2)];

                l=svgl(fine,p0,p1);

                pref=p1;

                R{id_line}=[R{id_line};l];
                
                num_row=num_row+1;
                
            elseif sta(end)=='h'
                p0=pref;

                p1=[str2double(P_line_inuse{num_row}),0]+pref;

                l=svgl(fine,p0,p1);

                pref=p1;

                R{id_line}=[R{id_line};l];
                
                num_row=num_row+1;
                
            elseif sta(end)=='C'
                p0=pref;

                p1=[str2double(P_line_inuse{num_row}),str2double(P_line_inuse{num_row+1})];

                p2=[str2double(P_line_inuse{num_row+2}),str2double(P_line_inuse{num_row+3})];

                p3=[str2double(P_line_inuse{num_row+4}),str2double(P_line_inuse{num_row+5})];

                pref=p3;

                b0123=svgc(fine,p0,p1,p2,p3);

                R{id_line}=[R{id_line};b0123];
                
                num_row=num_row+6;
                
                savep=p2;
                
            elseif sta(end)=='c'
                p0=pref;

                p1=[str2double(P_line_inuse{num_row}),str2double(P_line_inuse{num_row+1})]+pref;

                p2=[str2double(P_line_inuse{num_row+2}),str2double(P_line_inuse{num_row+3})]+pref;

                p3=[str2double(P_line_inuse{num_row+4}),str2double(P_line_inuse{num_row+5})]+pref;

                pref=p3;

                b0123=svgc(fine,p0,p1,p2,p3);

                R{id_line}=[R{id_line};b0123];
                
                num_row=num_row+6;
                
                savep=p2;

            elseif sta(end)=='S'
                p0=pref;

                if sta(end-1)=='C'||'c'||'S'||'s'
                    p1=2*pref-savep;
                else
                    p1=p0;
                end

                p2=[str2double(P_line_inuse{num_row}),str2double(P_line_inuse{num_row+1})];

                p3=[str2double(P_line_inuse{num_row+2}),str2double(P_line_inuse{num_row+3})];

                pref=p3;

                b0123=svgc(fine,p0,p1,p2,p3);

                R{id_line}=[R{id_line};b0123];
                
                num_row=num_row+4;
                
                savep=p2;
            
            elseif sta(end)=='s'
                p0=pref;

                if sta(end-1)=='C'||'c'||'S'||'s'
                    p1=2*pref-savep;
                else
                    p1=p0;
                end

                p2=[str2double(P_line_inuse{num_row}),str2double(P_line_inuse{num_row+1})]+pref;

                p3=[str2double(P_line_inuse{num_row+2}),str2double(P_line_inuse{num_row+3})]+pref;

                pref=p3;

                b0123=svgc(fine,p0,p1,p2,p3);

                R{id_line}=[R{id_line};b0123];
                
                num_row=num_row+4;
                
                savep=p2;
            
            end
        end
    end
    
%     R0=R{id_line};
%     R1=R0;
%     R2=[0,0;R0(1:end-1,:)];
% 
%     Rs=R1-R2;
% 
%     Rs=Rs(:,1) | Rs(:,2);
%     R{id_line}=R0(Rs==1,:);
end
%% plot svg

if graph
    fig=figure;

    for num_line=1:num_P_cell
        plot_R=R{num_line};
        outline=plot(plot_R(:,1),-plot_R(:,2));
        outline.Color=[0 0 0];
        outline.LineWidth=1.5;
        hold on
        axis equal
    end
end


end