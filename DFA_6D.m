% This matlab code undertakes the DF Analysis on a set of six categories
% and collapses the results onto a 3D point.
% Written by Martin Bencsik, Sep. 2015.
clear PCA_ranks
close all
highest_PCA_rank = 12;

PCA_ranks = 1:highest_PCA_rank;

% The following line requires two matlab m-files, 'dfa.m' and 'TWGen.m',
% supplied by Roy Goodacre from University of Manchester, UK.
[U,V, eigenval] = dfa([PCA_scores_01(PCA_ranks,:) PCA_scores_02(PCA_ranks,:) PCA_scores_03(PCA_ranks,:) PCA_scores_04(PCA_ranks,:) PCA_scores_05(PCA_ranks,:) PCA_scores_06(PCA_ranks,:)]', [zeros(1,size(PCA_scores_01,2)) ones(1,size(PCA_scores_02,2)) 2*ones(1,size(PCA_scores_03,2)) 3*ones(1,size(PCA_scores_04,2)) 4*ones(1,size(PCA_scores_05,2)) 5*ones(1,size(PCA_scores_06,2))],3);

DFA_spectrum_01 = sum(Eigenspectra(:,PCA_ranks)'.*(V(:,1)*ones(1,size(Eigenspectra,2))));
DFA_spectrum_02 = sum(Eigenspectra(:,PCA_ranks)'.*(V(:,2)*ones(1,size(Eigenspectra,2))));
DFA_spectrum_03 = sum(Eigenspectra(:,PCA_ranks)'.*(V(:,3)*ones(1,size(Eigenspectra,2))));

clear A_x A_y

A_x = sum((spectra_1).*(ones(size(spectra_1,1),1)*DFA_spectrum_01),2);
A_y = sum((spectra_1).*(ones(size(spectra_1,1),1)*DFA_spectrum_02),2);
A_z = sum((spectra_1).*(ones(size(spectra_1,1),1)*DFA_spectrum_03),2);

U_L = size(A_x,1);
centre_01 = [mean(A_x) mean(A_y) mean(A_z)];

A_x = [A_x; sum((spectra_2).*(ones(size(spectra_2,1),1)*DFA_spectrum_01),2)];
A_y = [A_y; sum((spectra_2).*(ones(size(spectra_2,1),1)*DFA_spectrum_02),2)];
A_z = [A_z; sum((spectra_2).*(ones(size(spectra_2,1),1)*DFA_spectrum_03),2)];

A_x = [A_x; sum((spectra_3).*(ones(size(spectra_3,1),1)*DFA_spectrum_01),2)];
A_y = [A_y; sum((spectra_3).*(ones(size(spectra_3,1),1)*DFA_spectrum_02),2)];
A_z = [A_z; sum((spectra_3).*(ones(size(spectra_3,1),1)*DFA_spectrum_03),2)];

A_x = [A_x; sum((spectra_4).*(ones(size(spectra_4,1),1)*DFA_spectrum_01),2)];
A_y = [A_y; sum((spectra_4).*(ones(size(spectra_4,1),1)*DFA_spectrum_02),2)];
A_z = [A_z; sum((spectra_4).*(ones(size(spectra_4,1),1)*DFA_spectrum_03),2)];

A_x = [A_x; sum((spectra_5).*(ones(size(spectra_5,1),1)*DFA_spectrum_01),2)];
A_y = [A_y; sum((spectra_5).*(ones(size(spectra_5,1),1)*DFA_spectrum_02),2)];
A_z = [A_z; sum((spectra_5).*(ones(size(spectra_5,1),1)*DFA_spectrum_03),2)];

A_x = [A_x; sum((spectra_6).*(ones(size(spectra_6,1),1)*DFA_spectrum_01),2)];
A_y = [A_y; sum((spectra_6).*(ones(size(spectra_6,1),1)*DFA_spectrum_02),2)];
A_z = [A_z; sum((spectra_6).*(ones(size(spectra_6,1),1)*DFA_spectrum_03),2)];

centre_02 = [mean(A_x((U_L+1):end)) mean(A_y((U_L+1):end)) mean(A_z((U_L+1):end))];


A = 1./sqrt((A_x - centre_01(1)).^2 + (A_y - centre_01(2)).^2 + (A_z - centre_01(3)).^2);


if (imag(A(1,1)) == 0)
    
    best_ranks = PCA_ranks
    figure(2)
    clf
    subplot(2,2,1)
    plot(A_x(1:size(spectra_1,1)),A_y(1:size(spectra_1,1)),'k.')
    hold on
    L_L = size(spectra_1,1) + 1;
    U_L = L_L + size(spectra_2,1)-1;
    plot(A_x(L_L:U_L),A_y(L_L:U_L),'r.');
    L_L = U_L + 1;
    U_L = L_L + size(spectra_3,1)-1;
    plot(A_x(L_L:U_L),A_y(L_L:U_L),'b.');
    L_L = U_L + 1;
    U_L = L_L + size(spectra_4,1)-1;
    plot(A_x(L_L:U_L),A_y(L_L:U_L),'g.');
    L_L = U_L + 1;
    U_L = L_L + size(spectra_5,1)-1;
    plot(A_x(L_L:U_L),A_y(L_L:U_L),'c.');
    L_L = U_L + 1;
    U_L = L_L + size(spectra_6,1)-1;
    plot(A_x(L_L:U_L),A_y(L_L:U_L),'m.');
    axis tight
    grid on
    xlabel('DF score 1')
    ylabel('DF score 2')
    
    subplot(2,2,2)
    plot(A_x(1:size(spectra_1,1)),A_z(1:size(spectra_1,1)),'k.')
    hold on
    L_L = size(spectra_1,1) + 1;
    U_L = L_L + size(spectra_2,1)-1;
    plot(A_x(L_L:U_L),A_z(L_L:U_L),'r.');
    L_L = U_L + 1;
    U_L = L_L + size(spectra_3,1)-1;
    plot(A_x(L_L:U_L),A_z(L_L:U_L),'b.');
    L_L = U_L + 1;
    U_L = L_L + size(spectra_4,1)-1;
    plot(A_x(L_L:U_L),A_z(L_L:U_L),'g.');
    L_L = U_L + 1;
    U_L = L_L + size(spectra_5,1)-1;
    plot(A_x(L_L:U_L),A_z(L_L:U_L),'c.');
    L_L = U_L + 1;
    U_L = L_L + size(spectra_6,1)-1;
    plot(A_x(L_L:U_L),A_z(L_L:U_L),'m.');
    axis tight
    grid on
    xlabel('DF score 1')
    ylabel('DF score 3')
    
    subplot(2,2,3)
    plot(A_y(1:size(spectra_1,1)),A_z(1:size(spectra_1,1)),'k.')
    hold on
    L_L = size(spectra_1,1) + 1;
    U_L = L_L + size(spectra_2,1)-1;
    plot(A_y(L_L:U_L),A_z(L_L:U_L),'r.');
    L_L = U_L + 1;
    U_L = L_L + size(spectra_3,1)-1;
    plot(A_y(L_L:U_L),A_z(L_L:U_L),'b.');
    L_L = U_L + 1;
    U_L = L_L + size(spectra_4,1)-1;
    plot(A_y(L_L:U_L),A_z(L_L:U_L),'g.');
    L_L = U_L + 1;
    U_L = L_L + size(spectra_5,1)-1;
    plot(A_y(L_L:U_L),A_z(L_L:U_L),'c.');
    L_L = U_L + 1;
    U_L = L_L + size(spectra_6,1)-1;
    plot(A_y(L_L:U_L),A_z(L_L:U_L),'m.');
    axis tight
    grid on
    xlabel('DF score 2')
    ylabel('DF score 3')
    
    subplot(2,2,4)
    
    plot(frequency_axis,DFA_spectrum_01,'linewidth',1)
    hold on
    plot(frequency_axis,DFA_spectrum_02,'r','linewidth',1)
    plot(frequency_axis,DFA_spectrum_03,'k','linewidth',1)
    title('Discriminant functions')
    ylabel('Acceleration')
    xlabel('Frequency (Hz)')
    axis tight
    grid on
    
end




