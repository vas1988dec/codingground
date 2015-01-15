function [ output_args ] = Untitled( input_args )
% in vecor se retin denumirile imaginilor pe care le vom deschide
% In imaginiInit se retin imaginile initiale de la 1 la 12
vector={'1.bmp';'2.bmp';'3.bmp';'4.bmp';'5.bmp'};
ImaginiInit=zeros(280,280,12);
for k=1:12
     nume=vector{k};
     Im=imread(nume);
ImaginiInit(:,:,k)=Im(:,:);
end;


% generez o matrice simetrica cu valori aleatoare, medie 0 si varianta
%4*4 - la fiecare executie este alta instanta de zgomot normrnd(mu,sigma,[m,n,...]) generates an m-by-n-by-... 
[n,m]=size(Im);
sigj=zeros(196,196);
for i=1:196
    sigj(i,i:196)=normrnd(0,4,[1,196-i+1]);
    for j=1:i-1
        sigj(i,j)=sigj(j,i);
    end;
end;

% Se vor liniariza imaginile initiale si se va calcula media pentru fiecare
% bloc liniarizat . In cazul nostru miu va avea 196 componte pentru fiecare
% bloc. Sigma(matricea de covarianta) este la fel
for i=1:n/14
    for j=1:n/14
        for k=1:196
            miu(k,i,j)=0;
        end;
        for k1=1:12
        for l=1:14
            for k=1:14
                T((l-1)*14+k,k1)=ImaginiInit((i-1)*14+l,(j-1)*14+k,k1);
                 T((l-1)*14+k,k1)=double( T((l-1)*14+k,k1))/255;
            end
        end
        %T(k1,1:196)=double(T(k1,1:196)/255);
        for k=1:196
        miu(k,i,j)=miu(k,i,j)+T(k,k1);
        end
        end
        miu(:,i,j)=miu(:,i,j)/12;
   
 sigma(:,:,i,j)=zeros(196,196);

       
        for k1=1:12
        for l=1:14
            for k=1:14
                 T((l-1)*14+k,k1)=ImaginiInit((i-1)*14+l,(j-1)*14+k,k1);
            end
        end
        sigma(:,:,i,j)=sigma(:,:,i,j)+(T(:,k1)-miu(:,i,j))*((T(:,k1)-miu(:,i,j)).');
        end

        sigma(:,:,i,j)=sigma(:,:,i,j)/(11);
    end
end
        
       
% matricea de covarianta a sgomotului este sigma sigj*sigj
   
% se vor perturba imaginile initiale liniarizate

for i=1:n/14
    for j=1:m/14
     
        for k1=1:12
        for l=1:14
            for k=1:14
                
                T(i,j,(l-1)*14+k,k1)=ImaginiInit((i-1)*14+l,(j-1)*14+k,k1);
            end
        end
        end
  
        vect=normrnd(0,1,[1,196]);
       
        vect=vect.';
       
        zgomot=sigj*vect;
         for k1=1:12
        for k=1:196
            PT(i,j,k,k1)=T(i,j,k,k1)+zgomot(k);
        end
         end
    end
end

% se vor reface imaginile initiale 
for i=1:n/14
    for j=1:m/14
        for k1=1:12
        for l=1:14
            for k=1:14
                PJ2((i-1)*14+l,(j-1)*14+k,k1)=PT(i,j,(l-1)*14+k,k1);
            end
        end
        end
    end
end

for k1=1:12
   B=uint8(PJ2(:,:,k1));
figure,imshow(B);  
end

%Imaginile perturbate se vor normaliza. De asemenea se va normaliza si
%sigmap
 for k1=1:12
ImgNorm(k1,1:m,1:n)=double(PJ2(1:m,1:n,k1))/255;
 end;
sigmap=sigj*sigj;
sigmap=sigmap/(255*255);

%se  liniarizeaza imaginile normalizate
for i=1:n/14
    for j=1:n/14
        for k1=1:12
        for t=1:14
            for p=1:14
                X(14*(t-1)+p,k1)=ImgNorm(k1,14*(i-1)+t,14*(j-1)+p);
            end
        end
%se vor centra valorile
        for v=1:196
        Y(v,k1)=X(v,k1)-miu(v,i,j);
        end
        end   
%lambda = eig(A,B) returns a vector containing the generalized eigenvalues of the pair, (A,B), that satisfy the equation, 

%Av = ?Bv

%where A and B are n-by-n matrices, v is a column vector of length n, and ? is a scalar. The values of ? that satisfy the equation are the generalized eigenvalues. 
%The corresponding values of v are generalized right eigenvectors.


%[V,D] = eig(___) returns two optional outputs for any of the previous syntaxes. V is a matrix whose columns are eigenvectors. D is a diagonal matrix containing the eigenvalues along the main diagonal. The form and scaling of V depends on the combination of input arguments:
%Matricea transformarii liniare
        [A,Lambda]=eig(sigmap,sigma(:,:,i,j));
        Z=zeros(196,12);
         for k1=1:12
% Aplica transformarea directa pentru decorelarea zgomotului
        Z(:,k1)=(A(:,:)')*Y(:,k1);
         end
% Aplica functia de contractie a codului pentru eliminarea zgomotului
         Z1=zeros(196,12);
         for k1=1:12
        for r=1:196
        if(Lambda(r,r)>0)
             Z1(r,k1)=sign(Z(r,k1))*max ([ 0 abs(Z(r,k1)) - sqrt(2)*((Lambda(r,r)))]);
        else 
            Z1(r,k1)=Z(r,k1);
        end
        end
         end
%Aplica transformarea inversa
         for k1=1:12
        Y0(:,k1)=pinv(A(:,:)')*Z1(:,k1);
         end
%Adaugarea mediei si recrearea imaginilor
         for k1=1:12
        for t=1:14
            for p=1:14
               R(14*(i-1)+t,14*(j-1)+p,k1)=Y0(14*(t-1)+p,k1)+miu(14*(t-1)+p,i,j);
            end
        end
         end
    end
end
% afisaarea imaginilor
for k1=1:12
B=uint8(R(:,:,k1)*255);
figure,imshow(B);  
end
end




