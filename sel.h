/*
 * void createLocalA(Matrix &A,mesh m){
    float u_bar = m.getParameter(ADJECTIVE_VELOCITY);
    A.at(0).at(0) += -u_bar/2;  A.at(0).at(1) += u_bar/2;
    A.at(1).at(0) += -u_bar/2;  A.at(1).at(1) += u_bar/2;
}

void createLocalB(Matrix &B,mesh m){
    float l = m.getParameter(ELEMENT_LENGTH);
    float nu = m.getParameter(DYNAMIC_VISCOSITY);
    B.at(0).at(0) += nu/l;      B.at(0).at(1) += -nu/l;
    B.at(1).at(0) += -nu/l;     B.at(1).at(1) += nu/l;
}

void createLocalC(Matrix &C,mesh m){
    float rho = m.getParameter(DENSITY);
    C.at(0).at(0) += -1/(2*rho);    C.at(0).at(1) += 1/(2*rho);
    C.at(1).at(0) += -1/(2*rho);    C.at(1).at(1) += 1/(2*rho);
}

void createLocalD(Matrix &D,mesh m){
    D.at(0).at(0) += -0.5;  D.at(0).at(1) += 0.5;
    D.at(1).at(0) += -0.5;  D.at(1).at(1) += 0.5;
}
*/
//MATRIX DE T
void createLocalA( Matrix &A, mesh m){
    float ct = m.getParameter(CT);
    A.at(0).at(0) += -ct/8;
    A.at(0).at(1) += ct/8;
    A.at(1).at(0) += -ct/8;
    A.at(1).at(1) += ct/8;
}
//MATRIZ DE K
void createLocalC(Matrix &C,mesh m){

    float l = m.getParameter(ELEMENT_LENGTH);
    float ck = m.getParameter(CK);
    C.at(0).at(0) += ck/l;
    C.at(0).at(1) += -ck/l;
    C.at(1).at(0) += -ck/l;
    C.at(1).at(1) += ck/l;
}
//MATRIZ DE LAMBDA
void createLocalD(Matrix &D,mesh m){
    float cl = m.getParameter(CL);
    D.at(0).at(0) += -cl/3;
    D.at(0).at(1) += cl/3;
    D.at(1).at(0) += -cl/3;
    D.at(1).at(1) += cl/3;
}
//MATRIZ DE V
void createLocalE(Matrix &E,mesh m){
    float l = m.getParameter(ELEMENT_LENGTH);
    float cv = m.getParameter(CV);
    E.at(0).at(0) += cv/l;
    E.at(0).at(1) += -cv/l;
    E.at(1).at(0) += -cv/l;
    E.at(1).at(1) += cv/l;
}
//MATRIZ ALFA
void createLocalF(Matrix &F,mesh m){

    float ca= m.getParameter(CA);
    F.at(0).at(0) += -3.0*ca/2;
    F.at(0).at(1) += 3.0*ca/2;
    F.at(1).at(0) += -3.0*ca/2;
    F.at(1).at(1) += 3.0*ca/2;

}
//MATRIZ DELTA
void createLocalG(Matrix &G,mesh m){
    float cd = m.getParameter(CD);
    G.at(0).at(0) += -cd/2.0;
    G.at(0).at(1) += cd/2.0;
    G.at(1).at(0) += -cd/2.0;
    G.at(1).at(1) += cd/2.0;
}

//a
Matrix createLocalK(int element,mesh &m){
    Matrix K,A,C,D,E,F,G;

    zeroes(A,2);
    zeroes(C,2);
    zeroes(D,2);
    zeroes(E,2);
    zeroes(F,2);
    zeroes(G,2);

    createLocalA(A,m);
    createLocalC(C,m);
    createLocalD(D,m);
    createLocalE(E,m);
    createLocalF(F,m);
    createLocalG(G,m);


    Vector row1, row2, row3, row4;

    row1.push_back(A.at(0).at(0)+C.at(0).at(0));
    row1.push_back(A.at(0).at(1)+C.at(0).at(1));
    row1.push_back(D.at(0).at(0)+E.at(0).at(0));
    row1.push_back(D.at(0).at(1)+E.at(0).at(1));

    row2.push_back(A.at(1).at(0)+C.at(1).at(0));
    row2.push_back(A.at(1).at(1)+C.at(1).at(1));
    row2.push_back(D.at(1).at(0)+E.at(1).at(0));
    row2.push_back(D.at(1).at(1)+E.at(1).at(1));

    row3.push_back(F.at(0).at(0));
    row3.push_back(F.at(0).at(1));
    row3.push_back(G.at(0).at(0));
    row3.push_back(G.at(0).at(1));

    row4.push_back(F.at(1).at(0));
    row4.push_back(F.at(1).at(1));
    row4.push_back(G.at(1).at(0));
    row4.push_back(G.at(1).at(1));

    K.push_back(row1);
    K.push_back(row2);
    K.push_back(row3);
    K.push_back(row4);
    K.push_back(row4);

    return K;
}

Vector createLocalb(int element,mesh &m){
    Vector b;

    float ctri = m.getParameter(CTRI);
    float l = m.getParameter(ELEMENT_LENGTH);
    float cn= m.getParameter(CN);

    b.push_back(ctri*l/2);
    b.push_back(ctri*l/2);
    b.push_back(cn*l/2);
    b.push_back(cn*l/2);

    return b;
}

void crearSistemasLocales(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        localKs.push_back(createLocalK(i,m));
        localbs.push_back(createLocalb(i,m));
    }
}

void assemblyK(element e,Matrix localK,Matrix &K,int nnodes){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = index1 + nnodes;
    int index4 = index2 + nnodes;

    K.at(index1).at(index1) += localK.at(0).at(0);
    K.at(index1).at(index2) += localK.at(0).at(1);
    K.at(index2).at(index1) += localK.at(1).at(0);
    K.at(index2).at(index2) += localK.at(1).at(1);

    K.at(index1).at(index3) += localK.at(0).at(2);
    K.at(index1).at(index4) += localK.at(0).at(3);
    K.at(index2).at(index3) += localK.at(1).at(2);
    K.at(index2).at(index4) += localK.at(1).at(3);

    K.at(index3).at(index1) += localK.at(2).at(0);
    K.at(index3).at(index2) += localK.at(2).at(1);
    K.at(index4).at(index1) += localK.at(3).at(0);
    K.at(index4).at(index2) += localK.at(3).at(1);

}

void assemblyb(element e,Vector localb,Vector &b){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;

    b.at(index1) += localb.at(0);
    b.at(index2) += localb.at(1);
}

void ensamblaje(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs,Matrix &K,Vector &b){
    int nnodes = m.getSize(NODES);
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        element e = m.getElement(i);
        assemblyK(e,localKs.at(i),K,nnodes);
        assemblyb(e,localbs.at(i),b);
    }
}


void applyDirichlet(mesh &m,Matrix &K,Vector &b){
    for(int i=0;i<m.getSize(DIRICHLET);i++){
        condition c = m.getCondition(i,DIRICHLET);

        int index = c.getNode1()-1;

        K.erase(K.begin()+index);
        b.erase(b.begin()+index);

        for(int row=0;row<K.size();row++){
            float cell = K.at(row).at(index);
            K.at(row).erase(K.at(row).begin()+index);
            b.at(row) += -1*c.getValue()*cell;
        }
    }
}


void calculate(Matrix &K, Vector &b, Vector &T){
    cout << "Iniciando calculo de respuesta...\n";
    Matrix Kinv;
    cout << "Calculo de inversa...\n";
    inverseMatrix(K,Kinv);
    cout << "Calculo de respuesta...\n";
    productMatrixVector(Kinv,b,T);
}