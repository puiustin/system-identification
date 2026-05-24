function DATA = make_IDSS()
%functia genereaza un obiect idss cu matrici aleatoare pentru testare

nx = 3;   % numarul de stari (ordinul sistemului)
nu = 5;   % numarul de semnale de intrare
ny = 3;   % numarul de semnale de iesire
Ts = 0.1; % perioada de esantionare (timpul de masura)

%generarea matricelor spatiului starilor cu valori random
A = randn(nx, nx); % matricea sistemului (evolutia starilor)
B = randn(nx, nu); % matricea de intrare (efectul intrarilor asupra starilor)
C = randn(ny, nx); % matricea de iesire (cum starile produc iesirile)
D = zeros(ny, nu); % matricea de legatura directa intre intrari si iesiri
K = randn(nx, ny); % matricea castigului kalman pentru procese de zgomot
X0 = zeros(nx,1);  % conditiile initiale ale starilor

%utilizarea constructorului idss pentru a crea obiectul modelului
DATA = idss(A,B,C,D,K,X0,Ts);

%configurarea numelui si notitelor pentru identificare
DATA.Name = 'IDSS_random';
DATA.Notes = 'model generat automat cu valori random pentru testare';
DATA.TimeUnit = 's'; % unitatea de timp este secunda
DATA.nk = zeros(1, nu); % intarzierea pe fiecare canal de intrare

%atribuirea de nume pentru fiecare stare (x1, x2, x3)
for i = 1:nx
    DATA.StateName{i} = ['x' num2str(i)];
end

%atribuirea de nume si unitati pentru semnalele de intrare (u1, u2...)
for i = 1:nu
    DATA.InputName{i} = ['u' num2str(i)];
    DATA.InputUnit{i} = 'unit_in';
end

%atribuirea de nume si unitati pentru semnalele de iesire (y1, y2...)
for i = 1:ny
    DATA.OutputName{i} = ['y' num2str(i)];
    DATA.OutputUnit{i} = 'unit_out';
end

%salvarea modelului creat intr-un fisier .mat pentru incarcare in ident
save('IDSS_random.mat','DATA');

%afisarea obiectului in consola pentru verificare
disp('model idss generat automat:');
disp(DATA);
end
