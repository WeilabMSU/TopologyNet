cd ./javaplex/;
load_javaplex;


LIGELE = {'C','N','O','S','CN','CO','CS','NO','NS','OS','CNO','CNS','COS','NOS','CNOS','CNOSPFClBrI','H','CH','NH','OH','SH','CNH','COH','CSH','NOH','NSH','OSH','CNOH','CNSH','COSH','NOSH','CNOSH','CNOSPFClBrIH','CCl','CClH','CBr','CBrH','CP','CF','CPH','CFH'};

formatSpec = '%f %f %f';
sizeA = [3,Inf];
formatSpecB = '%d %d';
sizeB = [2,Inf];


for j=1:41
    e2 = LIGELE{j};
    Name = strcat(pdb,'_',e2,'.pts');
    NameB = strcat(pdb,'_',e2,'.bds');
    OutName = strcat(pdb,'_',e2,'_level1.PH');
    if exist(strcat(DataDir,'/', Name), 'file') == 2
        fileID = fopen(strcat(DataDir,'/', Name), 'r');
        A = fscanf(fileID, formatSpec, sizeA);
        fclose(fileID);
        if exist(strcat(DataDir,'/', NameB), 'file') == 2
            fileID = fopen(strcat(DataDir,'/', NameB), 'r');
            B = fscanf(fileID, formatSpecB, sizeB);
            fclose(fileID);
        end
        distances = ones(size(A,2),size(A,2)).*100.0;
        for ii=1:size(A,2)
            for jj=(ii+1):size(A,2)
                dis = sqrt((A(1,ii) - A(1,jj))^2 + (A(2,ii) - A(2,jj))^2 + (A(3,ii) - A(3,jj))^2);
                distances(ii,jj) = dis;
                distances(jj,ii) = dis;
                distances(ii,ii) = 0.0;
            end
        end
        for ii=1:size(B,2)
            distances(B(1,ii),B(2,ii)) = 100.;
            distances(B(2,ii),B(1,ii)) = 100.;
        end
        m_space = metric.impl.ExplicitMetricSpace(distances);
        stream = api.Plex4.createVietorisRipsStream(m_space, 3, 10.0, 1000);
        persistence = api.Plex4.getModularSimplicialAlgorithm(3, 2);
        intervals = persistence.computeIntervals(stream);
        fileID = fopen(strcat(DataDir,'/',OutName), 'w');
        endpoints = homology.barcodes.BarcodeUtility.getEndpoints(intervals, 0, false);
        if size(endpoints,1) > 0
            %d = 0
            %size(endpoints,1)
            dims = zeros(1,size(endpoints,1));
            bars = [dims; endpoints(:,1)'; endpoints(:,2)'];
            fprintf(fileID, '%d %4.4f %4.4f\n', bars);
        end
        endpoints = homology.barcodes.BarcodeUtility.getEndpoints(intervals, 1, false);
        if size(endpoints,1) > 0
            %d = 1
            %size(endpoints,1)
            dims = ones(1,size(endpoints,1));
            bars = [dims; endpoints(:,1)'; endpoints(:,2)'];
            fprintf(fileID, '%d %4.4f %4.4f\n', bars);
        end
        endpoints = homology.barcodes.BarcodeUtility.getEndpoints(intervals, 2, false);
        if size(endpoints,1) > 0
            %d = 2
            %size(endpoints,1)
            dims = ones(1,size(endpoints,1))*2;
            bars = [dims; endpoints(:,1)'; endpoints(:,2)'];
            fprintf(fileID, '%d %4.4f %4.4f\n', bars);
        end
        fclose(fileID);
    end
end

exit
