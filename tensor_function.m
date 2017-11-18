function tensor_function(Index)

radius = 1;
pradius = 1;
ppatchsize = pradius*2+1;
mradius = radius + pradius + 1; % to ensure patches extracted do not go out of image region
patchsize = radius*2+1;

Base = [];
for i=(-radius):radius
    for j = (-radius):radius
        for k = (-radius):radius
            if i==0&&j==0&&k==0
                continue;
            end
            Base = [Base [i;j;k]];
        end
    end
end
dist = sqrt(sum(Base.*Base));
Base = Base./repmat(dist,3,1);

Base6 = cell(patchsize^3-1);
for i=1:(patchsize^3-1)
    Base6{i} = Base(:,i)*(Base(:,i)');
end


tic

baseFolder=['/Volumes/neural-network-registration-data/NC/FunImgARCF/',Index,'/'];
mrFilename='Filtered_4DVolume.nii';

[mrData,mrHead]=y_Read([baseFolder,mrFilename]);
mrData(isnan(mrData)==1)=0;mrData=squeeze(mrData);
[maskData,~]=y_Read(['/Volumes/neural-network-registration-data/NC/T1Img_ori_seg/',Index,'/rc1T1_crop.nii']);

[rows,cols,depths,frames]=size(mrData);
eigenVecMat1=zeros(rows,cols,depths,3);
eigenVecMat2=zeros(rows,cols,depths,3);
eigenVecMat3=zeros(rows,cols,depths,3);
tensorMat=zeros(rows,cols,depths,6);
eigenValueMat=zeros(rows,cols,depths,3);
faMat=zeros(rows,cols,depths,1);
rdMat=zeros(rows,cols,depths,1);
mdMat=zeros(rows,cols,depths,1);
vrMat=zeros(rows,cols,depths,1);
% v = (-radius):radius;
pv = (-pradius):(pradius);

for row=mradius:rows-mradius
    for col=mradius:cols-mradius
        for depth=mradius:depths-mradius
            if maskData(row,col,depth)>0%do tensor
                B=zeros(frames,ppatchsize^3,patchsize^3);

                for i = 1:patchsize
                    for j = 1:patchsize
                        for k = 1:patchsize
                            Tmp=permute(mrData(row-radius-1+i+pv,col-radius-1+j+pv,depth-radius-1+k+pv,:),[4 3 2 1]);
                            B(:,:,(i-1)*patchsize*patchsize+(j-1)*patchsize+k) = squeeze(reshape(Tmp,[],frames,ppatchsize^3));
                        end
                    end
                end
%                A = mrData(row+v,col+v,depth+v,:); % 3x3x3xtime

%                 Tmp = permute(A,[4 3 2 1]);
%                 B=squeeze(reshape(Tmp,[],size(A,4), patchsize^3)); %
                MM=zeros(patchsize^3-1,1);
                for i=1:ppatchsize^3
                    MM = MM+abs(corr(squeeze(B(:,i,(patchsize^3+1)/2)),squeeze(B(:,i,[1:(patchsize^3-1)/2 (patchsize^3+3)/2:end]))))';
                end
                MM=MM/ppatchsize^3;
                MM(isnan(MM)) = 0;
                
                currMat = zeros(3,3);
                for i=1:(patchsize^3-1)
                    currMat = currMat + MM(i)*Base6{i}*MM(i);
                end
                tensorMat(row,col,depth,:)=[currMat(1),currMat(2),currMat(3),currMat(5),currMat(6),currMat(9)];
                [coeff,latent]=pcacov(currMat);%it is equal to get eigen value for currMat'*currMat
                if size(coeff,2)<3
%                     fprintf('size of coeff is not 3\n');
                    continue;
                end
                eigenVecMat1(row,col,depth,:)=coeff(:,1);
                eigenVecMat2(row,col,depth,:)=coeff(:,2);
                eigenVecMat3(row,col,depth,:)=coeff(:,3);
                eigenValueMat(row,col,depth,:)=latent;
                meanEigen=mean(latent);
                
                fa=sqrt(2/3*((latent(1)-meanEigen)^2+(latent(2)-meanEigen)^2+(latent(3)-meanEigen)^2)/(latent(1)^2+latent(2)^2+latent(3)^2));
                md=meanEigen;
                rd=(latent(2)+latent(3))/2;
                vr=latent(1)*latent(2)*latent(3)*sqrt(meanEigen*meanEigen*meanEigen);
                faMat(row,col,depth)=fa;
                mdMat(row,col,depth)=md;
                rdMat(row,col,depth)=rd;
                vrMat(row,col,depth)=vr;
 
            end
        end
    end
end
toc
rest_Write4DNIfTI(tensorMat,mrHead,[baseFolder,'tensor.nii']);
rest_WriteNiftiImage(faMat(:,:,:),mrHead,[baseFolder,'fa.nii']);
rest_WriteNiftiImage(mdMat(:,:,:),mrHead,[baseFolder,'md.nii']);
rest_WriteNiftiImage(rdMat(:,:,:),mrHead,[baseFolder,'rd.nii']);
rest_WriteNiftiImage(vrMat(:,:,:),mrHead,[baseFolder,'vr.nii']);
rest_Write4DNIfTI(eigenVecMat1,mrHead,[baseFolder,'eigenMat_1.nii']);
rest_Write4DNIfTI(eigenVecMat2,mrHead,[baseFolder,'eigenMat_2.nii']);
rest_Write4DNIfTI(eigenVecMat3,mrHead,[baseFolder,'eigenMat_3.nii']);
rest_Write4DNIfTI(eigenValueMat,mrHead,[baseFolder,'eigenValue.nii']);
quit
