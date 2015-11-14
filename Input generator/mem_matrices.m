function [Amem,Bmem,Cmem,sizechi,dof_mem] = mem_matrices(cpl_mem,Ar,Br,Cr)
% [Amem,Bmem,Cmem,sizechi,dof_mem] = mem_matrices(cpl_mem,Ar,Br,Cr)


warning('off')
% Calculation of number of degrees of freedom, and of how many times each of them is considered
dof_mem = zeros(6,2);
for k1 = 1:size(cpl_mem,1)
    cdof = cpl_mem(k1,1);
    if dof_mem(cdof,1) == 0
        dof_mem(cdof,1) = 1;
        dof_mem(cdof,2) = 1;
    else
        dof_mem(cdof,2) = dof_mem(cdof,2) + 1;
    end
end

if nargin == 1
    % State space matrices loading
    sizechi = 0;
    Amem = [];
    Bmem = [];
    Cmem = [];
    for k1 = 1:size(cpl_mem,1)
        eval([' load Ar' num2str(cpl_mem(k1,1)) num2str(cpl_mem(k1,2)) '.txt -ascii; '])
        eval([' load Br' num2str(cpl_mem(k1,1)) num2str(cpl_mem(k1,2)) '.txt -ascii; '])
        eval([' load Cr' num2str(cpl_mem(k1,1)) num2str(cpl_mem(k1,2)) '.txt -ascii; '])
        eval([' nA = Ar' num2str(cpl_mem(k1,1)) num2str(cpl_mem(k1,2)) '; '])
        eval([' nB = Br' num2str(cpl_mem(k1,1)) num2str(cpl_mem(k1,2)) '; '])
        eval([' nC = Cr' num2str(cpl_mem(k1,1)) num2str(cpl_mem(k1,2)) '; '])
        nlin1 = size(Amem,1);
        ncol1 = size(Amem,2);
        nlin2 = size(nA,1);
        ncol2 = size(nA,2);
        Amem = [Amem zeros(nlin1,ncol2);zeros(nlin2,ncol1) nA];
        if k1 == 1
            Bmem = nB;
            Cmem = nC;
        else
            nlinb1 = size(Bmem,1);
            ncolb1 = size(Bmem,2);
            nlinb2 = size(nB,1);
            ncolb2 = size(nB,2);
            nlinc1 = size(Cmem,1);
            ncolc1 = size(Cmem,2);
            nlinc2 = size(nC,1);
            ncolc2 = size(nC,2);
            if cpl_mem(k1,1) == cpl_mem(k1-1,1)
                Bmem = [Bmem;zeros(nlinb2,ncolb1-1) nB];
                Cmem = [Cmem [zeros(nlinc1-1,ncolc2);nC]];
            else
                Bmem = [Bmem zeros(nlinb1,1);zeros(nlinb2,ncolb1) nB];
                Cmem = [Cmem zeros(nlinc1,ncolc2);zeros(1,ncolc1) nC];
            end
        end
    end
    sizechi = size(Amem,1);
else
    
    sizechi = 0;
    Amem = [];
    Bmem = [];
    Cmem = [];
    for k1 = 1:size(cpl_mem,1)
        nA = Ar{cpl_mem(k1,1),cpl_mem(k1,2)};
        nB = Br{cpl_mem(k1,1),cpl_mem(k1,2)};
        nC = Cr{cpl_mem(k1,1),cpl_mem(k1,2)};
        nlin1 = size(Amem,1);
        ncol1 = size(Amem,2);
        nlin2 = size(nA,1);
        ncol2 = size(nA,2);
        Amem = [Amem zeros(nlin1,ncol2);zeros(nlin2,ncol1) nA];
        if k1 == 1
            Bmem = nB;
            Cmem = nC;
        else
            nlinb1 = size(Bmem,1);
            ncolb1 = size(Bmem,2);
            nlinb2 = size(nB,1);
            ncolb2 = size(nB,2);
            nlinc1 = size(Cmem,1);
            ncolc1 = size(Cmem,2);
            nlinc2 = size(nC,1);
            ncolc2 = size(nC,2);
            if cpl_mem(k1,1) == cpl_mem(k1-1,1)
                Bmem = [Bmem;zeros(nlinb2,ncolb1-1) nB];
                Cmem = [Cmem [zeros(nlinc1-1,ncolc2);nC]];
            else
                Bmem = [Bmem zeros(nlinb1,1);zeros(nlinb2,ncolb1) nB];
                Cmem = [Cmem zeros(nlinc1,ncolc2);zeros(1,ncolc1) nC];
            end
        end
    end
    sizechi = size(Amem,1);
end
warning('on')