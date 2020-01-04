% %graphs for Exercise 3 PDS
% datav0 = importdata('datav0.csv');
% datav1 = importdata('datav1.csv');
% datav2 = importdata('datav2.csv');
% datav3 = importdata('datav3.csv');
% 
% %graphs for constant n and k varius
% %n=100
%     v0_n100= datav0(57:62,:);
%     v1_n100= datav1(57:62,:);
%     v2_n100= datav2(57:62,:);
%     v3_n100= datav3(57:62,:);
%     figure(1);
%     plot(v0_n100(:,2),v0_n100(:,3),'b');
%     hold on;
%     plot(v0_n100(:,2),v1_n100(:,3),'g');
%     hold on;
%     plot(v0_n100(:,2),v2_n100(:,3),'r');
%     hold on;
%     plot(v0_n100(:,2),v3_n100(:,3),'magenta');
%     hold on;
%     title("Execution times for n=100 and various k.(V0-V3)");
%     xlabel("k: steps");
%     ylabel("Execution time (seconds)");
% 
%     legend("V0", "V1","V2","V3 ");
%     
%     figure(2);
%     plot(v0_n100(:,2),v1_n100(:,3),'g');
%     hold on;
%     plot(v0_n100(:,2),v2_n100(:,3),'r');
%     hold on;
%     plot(v0_n100(:,2),v3_n100(:,3),'magenta');
%     hold on;
%     title("Execution times for n=100 and various k.(V1-V3)");
%     xlabel("k: steps");
%     ylabel("Execution time (seconds)");
% 
%     legend("V1","V2","V3 ");
% 
% 
% % %n=250
% % 
% % %n=500
%     v0_n500= datav0(43:49,:);
%     v1_n500= datav1(43:49,:);
%     v2_n500= datav2(43:49,:);
%     v3_n500= datav3(43:49,:);
%     figure(3);
%     plot(v0_n500(:,2),v0_n500(:,3),'b');
%     hold on;
%     plot(v0_n500(:,2),v1_n500(:,3),'g');
%     hold on;
%     plot(v0_n500(:,2),v2_n500(:,3),'r');
%     hold on;
%     plot(v0_n500(:,2),v3_n500(:,3),'magenta');
%     hold on;
%     title("Execution times for n=500 and various k.(V0-V3)");
%     xlabel("k: steps");
%     ylabel("Execution time (seconds)");
% 
%     legend("V0", "V1","V2","V3 ");
%     
%     figure(4);
%     plot(v0_n500(:,2),v1_n500(:,3),'g');
%     hold on;
%     plot(v0_n500(:,2),v2_n500(:,3),'r');
%     hold on;
%     plot(v0_n500(:,2),v3_n500(:,3),'magenta');
%     hold on;
%     title("Execution times for n=500 and various k.(V1-V3)");
%     xlabel("k: steps");
%     ylabel("Execution time (seconds)");
% 
%     legend("V1","V2","V3 ");
% 
% 
% % %n=1000
%     v0_n1000= datav0(29:35,:);
%     v1_n1000= datav1(29:35,:);
%     v2_n1000= datav2(29:35,:);
%     v3_n1000= datav3(29:35,:);
%     figure(5);
%     plot(v0_n1000(:,2),v0_n1000(:,3),'b');
%     hold on;
%     plot(v0_n1000(:,2),v1_n1000(:,3),'g');
%     hold on;
%     plot(v0_n1000(:,2),v2_n1000(:,3),'r');
%     hold on;
%     plot(v0_n1000(:,2),v3_n1000(:,3),'magenta');
%     hold on;
%     title("Execution times for n=1000 and various k.(V0-V3)");
%     xlabel("k: steps");
%     ylabel("Execution time (seconds)");
% 
%     legend("V0", "V1","V2","V3 ");
%     
%     figure(6);
%     plot(v0_n1000(:,2),v1_n1000(:,3),'g');
%     hold on;
%     plot(v0_n1000(:,2),v2_n1000(:,3),'r');
%     hold on;
%     plot(v0_n1000(:,2),v3_n1000(:,3),'magenta');
%     hold on;
%     title("Execution times for n=1000 and various k.(V1-V3)");
%     xlabel("k: steps");
%     ylabel("Execution time (seconds)");
% 
%     legend("V1","V2","V3 ");
% 
% %n=5000
%     v0_n5000= datav0(1:7,:);
%     v1_n5000= datav1(1:7,:);
%     v2_n5000= datav2(1:7,:);
%     v3_n5000= datav3(1:7,:);
%     
%     figure(7);
%     plot(v0_n5000(:,2),v0_n5000(:,3),'b');
%     hold on;   
%     plot(v0_n5000(:,2),v1_n5000(:,3),'g');
%     hold on;
%     plot(v0_n5000(:,2),v2_n5000(:,3),'r');
%     hold on;
%     plot(v0_n5000(:,2),v3_n5000(:,3),'magenta');
%     hold on;
%     title("Execution times for n=5000 and various k.(V0-V3)");
%     xlabel("k: steps");
%     ylabel("Execution time (seconds)");
% 
%     legend("V0", "V1","V2","V3 ");
%     
%     figure(8);
%     plot(v0_n5000(:,2),v1_n5000(:,3),'g');
%     hold on;
%     plot(v0_n5000(:,2),v2_n5000(:,3),'r');
%     hold on;
%     plot(v0_n5000(:,2),v3_n5000(:,3),'magenta');
%     hold on;
%     title("Execution times for n=5000 and various k.(V1-V3)");
%     xlabel("k: steps");
%     ylabel("Execution time (seconds)");
% 
%     legend( "V1","V2","V3 ");


%graphs for constant k and n varius
%First I sort the datav0-datav3 matrices per 2nd comlumn
    %κ=1
    v0_k1= datav0(1:8,:);
    v1_k1= datav1(1:8,:);
    v2_k1= datav2(1:8,:);
    v3_k1= datav3(1:8,:);
    
    figure(9);
    plot(v0_k1(:,1),v0_k1(:,3),'b');
    hold on;   
    plot(v0_k1(:,1),v1_k1(:,3),'g');
    hold on;
    plot(v0_k1(:,1),v2_k1(:,3),'r');
    hold on;
    plot(v0_k1(:,1),v3_k1(:,3),'magenta');
    hold on;
    title("Execution times for k=1.(V0-V3)");
    xlabel("n: one dimension of the nxn Lattice");
    ylabel("Execution time (seconds)");

    legend("V0", "V1","V2","V3 ");
    
    figure(10);
    plot(v0_k1(:,1),v1_k1(:,3),'g');
    hold on;
    plot(v0_k1(:,1),v2_k1(:,3),'r');
    hold on;
    plot(v0_k1(:,1),v3_k1(:,3),'magenta');
    hold on;
    title("Execution times for k=1.(V1-V3)");
    xlabel("n: one dimension of the nxn Lattice");
    ylabel("Execution time (seconds)");

    legend( "V1","V2","V3 ");
    
    
    %κ=10
    v0_k1= datav0(19:27,:);
    v1_k1= datav1(19:27,:);
    v2_k1= datav2(19:27,:);
    v3_k1= datav3(19:27,:);
    
    figure(11);
    plot(v0_k1(:,1),v0_k1(:,3),'b');
    hold on;   
    plot(v0_k1(:,1),v1_k1(:,3),'g');
    hold on;
    plot(v0_k1(:,1),v2_k1(:,3),'r');
    hold on;
    plot(v0_k1(:,1),v3_k1(:,3),'magenta');
    hold on;
    title("Execution times for k=10.(V0-V3)");
    xlabel("n: one dimension of the nxn Lattice");
    ylabel("Execution time (seconds)");

    legend("V0", "V1","V2","V3 ");
    
    figure(12);
    plot(v0_k1(:,1),v1_k1(:,3),'g');
    hold on;
    plot(v0_k1(:,1),v2_k1(:,3),'r');
    hold on;
    plot(v0_k1(:,1),v3_k1(:,3),'magenta');
    hold on;
    title("Execution times for k=10(V1-V3)");
    xlabel("n: one dimension of the nxn Lattice");
    ylabel("Execution time (seconds)");

    legend( "V1","V2","V3 ");
    
    %κ=50
    v0_k1= datav0(46:54,:);
    v1_k1= datav1(46:54,:);
    v2_k1= datav2(46:54,:);
    v3_k1= datav3(46:54,:);
    
    figure(13);
    plot(v0_k1(:,1),v0_k1(:,3),'b');
    hold on;   
    plot(v0_k1(:,1),v1_k1(:,3),'g');
    hold on;
    plot(v0_k1(:,1),v2_k1(:,3),'r');
    hold on;
    plot(v0_k1(:,1),v3_k1(:,3),'magenta');
    hold on;
    title("Execution times for k=50.(V0-V3)");
    xlabel("n: one dimension of the nxn Lattice");
    ylabel("Execution time (seconds)");

    legend("V0", "V1","V2","V3 ");
    
    figure(14);
    plot(v0_k1(:,1),v1_k1(:,3),'g');
    hold on;
    plot(v0_k1(:,1),v2_k1(:,3),'r');
    hold on;
    plot(v0_k1(:,1),v3_k1(:,3),'magenta');
    hold on;
    title("Execution times for k=50(V1-V3)");
    xlabel("n: one dimension of the nxn Lattice");
    ylabel("Execution time (seconds)");

    legend( "V1","V2","V3 ");
    
    %κ=100
    v0_k1= datav0(55:63,:);
    v1_k1= datav1(55:63,:);
    v2_k1= datav2(55:63,:);
    v3_k1= datav3(55:63,:);
    
    figure(15);
    plot(v0_k1(:,1),v0_k1(:,3),'b');
    hold on;   
    plot(v0_k1(:,1),v1_k1(:,3),'g');
    hold on;
    plot(v0_k1(:,1),v2_k1(:,3),'r');
    hold on;
    plot(v0_k1(:,1),v3_k1(:,3),'magenta');
    hold on;
    title("Execution times for k=100.(V0-V3)");
    xlabel("n: one dimension of the nxn Lattice");
    ylabel("Execution time (seconds)");

    legend("V0", "V1","V2","V3 ");
    
    figure(16);
    plot(v0_k1(:,1),v1_k1(:,3),'g');
    hold on;
    plot(v0_k1(:,1),v2_k1(:,3),'r');
    hold on;
    plot(v0_k1(:,1),v3_k1(:,3),'magenta');
    hold on;
    title("Execution times for k=100(V1-V3)");
    xlabel("n: one dimension of the nxn Lattice");
    ylabel("Execution time (seconds)");

    legend( "V1","V2","V3 ");



%3D plot 
    % scatter3(datav0(:,1),datav0(:,2),datav0(:,3),'filled');
    % hold on;
    % scatter3(datav1(:,1),datav1(:,2),datav1(:,3),'filled');
    % hold on;
    % scatter3(datav2(:,1),datav2(:,2),datav2(:,3),'filled');
    % hold on;
    % scatter3(datav3(:,1),datav3(:,2),datav3(:,3),'filled');
