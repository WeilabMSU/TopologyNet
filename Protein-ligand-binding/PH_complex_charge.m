cd ./javaplex/;
load_javaplex;

PROELE = {'C','N','O','S','H'};
LIGELE = {'C','N','O','S','P','F','Cl','Br','I','H'};

formatSpec = '%d %f %f %f %f';
sizeA = [5,Inf];


for j=1:5
    for k=1:10
        e1 = PROELE{j}; e2 = LIGELE{k};
        Name = strcat(pdb,'_',e1,'_',e2,'_16.0_chg.pts');
        OutName = strcat(pdb,'_',e1,'_',e2,'_16.0_chg.PH');
        if exist(strcat(DataDir,'/', Name), 'file') == 2
            fileID = fopen(strcat(DataDir,'/', Name), 'r');
            A = fscanf(fileID, formatSpec, sizeA);
            fclose(fileID);
            distances = ones(size(A,2),size(A,2)).*100.0;
            for ii=1:size(A,2)
                for jj=(ii+1):size(A,2)
                    if A(1,ii)+A(1,jj) == 1
                        dis = sqrt((A(2,ii) - A(2,jj))^2 + (A(3,ii) - A(3,jj))^2 + (A(4,ii) - A(4,jj))^2);
                        if dis < 0.002
                            dis
                        end
                        qi = A(5,ii);
                        qj = A(5,jj);
                        distances(ii,jj) = 1.0/(1.0+exp(-100.0*qi*qj/dis));
                        distances(jj,ii) = 1.0/(1.0+exp(-100.0*qi*qj/dis));
                        %distances(ii,jj) = dis;
                        %distances(jj,ii) = dis;
                        distances(ii,ii) = 0.0;
                    end
                end
            end
            m_space = metric.impl.ExplicitMetricSpace(distances);
            stream = api.Plex4.createVietorisRipsStream(m_space, 1, 2.0, 20000);
            persistence = api.Plex4.getModularSimplicialAlgorithm(1, 2);
            intervals = persistence.computeIntervals(stream);
            endpoints = homology.barcodes.BarcodeUtility.getEndpoints(intervals, 0, false);
            dims = zeros(1,size(endpoints,1));
            bars = [dims; endpoints(:,1)'; endpoints(:,2)'];
            fileID = fopen(strcat(DataDir,'/', OutName), 'w');
            fprintf(fileID, '%d %4.4f %4.4f\n', bars);
            fclose(fileID);
            clear m_space; clear stream; clear persistence; clear intervals; clear endpoints;
        end
    end
end

exit
