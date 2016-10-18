# funicular

clear
N       = 35;     // Nombre de ville
bet     = 1000;      // Inverse température dans Métropolis
maxiter = 1000;    // Maximum d'itération dans Métroplis


fileVdC = 'VdC'+string(N)+'.txt";
[stat]=unix("ls "+fileVdC);

if stat~=0 then
    points = grand(N,2,"unf",0,1);
    fprintfMat("VdC'+string(N)+'.txt",points);
else
    points = fscanfMat(fileVdC);
end

//points=grand(N,2,"unf",0,1);
//fprintfMat("VdC'+string(N)+'.txt",points)

points = fscanfMat("VdC'+string(N)+'.txt");

scf(0); clf; plot(points(:,1),points(:,2),'r*')
a=gca(); a.isoview='on';

// Calcule la longueur du chemin indexé par la permutation ind
function [longueur]=H(pts)
    ptsc = [pts ; pts(1,:)];
    longueur = sqrt(sum(diff(ptsc,1,1).^2));
endfunction

// trace le chemin ind
function []=tracer(pts,ind)
    plot(pts(ind,1),pts(ind,2))
endfunction

// 1ere matrice de sélection : on échange 2 points au hasard sur le chemin
function [Xk1,HX1] = simulP2(Xk)
    // 1ère étape : on tire Y selon P
    Yk1 = Xk;
    kk = grand(2,1,'uin',1,size(Xk,1));
    while diff(kk)==0,
        kk = grand(2,1,'uin',1,size(Xk,1));
    end
    Yk1(kk(1))=Xk(kk(2));
    Yk1(kk(2))=Xk(kk(1));

    Hy = H(points(Yk1,:)); Hx = H(points(Xk,:));
    // 2ème étape : on accepte/rejette
    if Hy<=Hx then
        Xk1 = Yk1;
        HX1 = Hy;
    else
        a = grand(1,1,'unf',0,1);
        if a<exp(-bet*(Hy-Hx)) then
            Xk1 = Yk1;
            HX1 = Hy;
        else
            Xk1 = Xk;
            HX1 = Hx;
        end
    end
endfunction

// 2ere matrice de sélection : on échange 2 points au hasard sur le chemin
// et on inverse le chemin entre les deux
function [Xk1,HX1] = simulP4(Xk)
    // 1ère étape : on tire Y selon P
    Yk1 = Xk;
    kk = grand(2,1,'uin',1,size(Xk,1));
    while diff(kk)==0,
        kk = grand(2,1,'uin',1,size(Xk,1));
    end
    Yk1(kk(1))=Xk(kk(2));
    Yk1(kk(2))=Xk(kk(1));

    Yk1(min(kk)+1:max(kk)-1)=Xk(max(kk)-1:-1:min(kk)+1);

    Hy = H(points(Yk1,:)); Hx = H(points(Xk,:));
    // 2ème étape : on accepte/rejette
    if Hy<=Hx then
        Xk1 = Yk1;
        HX1 = Hy;
    else
        a = grand(1,1,'unf',0,1);
        if a<exp(-bet*(Hy-Hx)) then
            Xk1 = Yk1;
            HX1 = Hy;
        else
            Xk1 = Xk;
            HX1 = Hx;
        end
    end
endfunction

/////////////////////////
//  ALGO DE METROPOLIS
/////////////////////////
// Tirage initiale (loi uniforme sur le groupe des permutations)H
[vv,X0] = gsort(grand(N,1,'unf',0,1),'lr','i');
HX0 = H(points(X0,:));

iter=1;
tab(1)=[HX0];
while iter<=maxiter,

    [X1,HX1] = simulP4(X0);

    tab(iter + 1) = sqrt(HX1);
    iter = iter + 1;
    X0 = X1;
end

scf(0); clf;
xtitle("Recuit simulé, beta = "+string(bet)+", N = "+string(N)+", MaxIter = "+string(maxiter)+", H = "+string(sqrt(H(points(X0,:)))));
//subplot(211)
plot(tab);
//xtitle("Longueur^2 des chemins","Iterations","H(X)")

X0c = [X0;X0(1)];
//subplot(212) 
scf(1); clf;
plot(points(X0c,1),points(X0c,2))
plot(points(:,1),points(:,2),'r*')
a=gca(); a.isoview='on';
//a.data_bounds=[0,0;1,1];
//xtitle("Dernier cheminement obtenu -- Iter = "+string(maxiter))
xtitle("Recuit simulé, beta = "+string(bet)+", N = "+string(N)+", MaxIter = "+string(maxiter)+", H = "+string(sqrt(H(points(X0,:)))));

fprintfMat("VdC"+string(N)+"_bet"+string(bet)+"_MaxIter"+string(maxiter)+".txt",X0);



