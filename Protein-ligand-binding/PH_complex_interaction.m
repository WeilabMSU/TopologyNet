cd ./javaplex/;
load_javaplex;

PROELE = {'C','N','O','S'};
LIGELE = {'C','N','O','S','P','F','Cl','Br','I'};

formatSpec = '%d %f %f %f';
sizeA = [4,Inf];

for j=1:4
    for k=1:9
        e1 = PROELE{j}; e2 = LIGELE{k};
        Name = strcat(pdb,'_',e1,'_',e2,'_50.0.pts');
        OutName = strcat(pdb,'_',e1,'_',e2,'_50.0_interaction.PH');
        if exist(strcat(DataDir, '/', Name), 'file') == 2
            fileID = fopen(strcat(DataDir, '/', Name), 'r');
            A = fscanf(fileID, formatSpec, sizeA);
            fclose(fileID);
            distances = ones(size(A,2),size(A,2)).*100.0;
            for ii=1:size(A,2)
                for jj=(ii+1):size(A,2)
                    if A(1,ii)+A(1,jj) == 1
                        dis = sqrt((A(2,ii) - A(2,jj))^2 + (A(3,ii) - A(3,jj))^2 + (A(4,ii) - A(4,jj))^2);
                        distances(ii,jj) = dis;
                        distances(jj,ii) = dis;
                        distances(ii,ii) = 0.0;
                    end
                end
            end
            m_space = metric.impl.ExplicitMetricSpace(distances);
            stream = api.Plex4.createVietorisRipsStream(m_space, 1, 50.0, 1000);
            persistence = api.Plex4.getModularSimplicialAlgorithm(1, 2);
            intervals = persistence.computeIntervals(stream);
            endpoints = homology.barcodes.BarcodeUtility.getEndpoints(intervals, 0, false);
            dims = zeros(1,size(endpoints,1));
            bars = [dims; endpoints(:,1)'; endpoints(:,2)'];
            fileID = fopen(strcat(DataDir, '/', OutName), 'w');
            fprintf(fileID, '%d %4.4f %4.4f\n', bars);
            fclose(fileID);
            clear m_space;
            clear stream;
            clear persistence;
            clear intervals;
            clear endpoints;
        end
    end
end

exit
