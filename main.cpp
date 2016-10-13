#include "cplex.h"

void timeused(double *time)
{
	static double tstart, tend, tprev;

	if (time == NULL) {
		clock(); /* one extra call to initialize clock */
		tstart = tprev = clock();
	} else {
		tend = clock();
		if (tend < tprev) 
            tstart -= ULONG_MAX; /* wraparound occured */
		tprev = tend;
		*time = (tend-tstart) / CLOCKS_PER_SEC; /* convert to seconds */
	}
}
 

int main(int argc, char *argv[]){

	char *inst; 
	if(argc < 2)
		inst = "problema.txt";
	else
		inst = argv[1];

		
	Problem *modelo = new Problem(inst); //pelo cmd
	DFS<Problem> search( modelo );
	double time1, time2;
    ofstream arq;
    arq.open("padroes.txt");

    int contador = 0;
    Problem *sol;
    list< list<float> > padrao;

    int c;
    bool bo = false;

    timeused(NULL);
    while ((sol = search.next())) {
        //impressao do padrao vazio com o #tamanhos do tipo 0
        if(!bo){
            c = sol->get_k0();
            bo = true;
        }
        list<float> aux;
        //sol->imprimir();
        sol->imprimir_lista(aux);
        padrao.push_back(aux);
        ++ contador;
        delete sol;
        
    }
	//contabilizar tempo gasto na geracao de padroes
    timeused(&time1);

	
	//numero de padroes
    arq << contador+1 << endl;
	
	//imprimir padrao 0
    arq << 0 << " ";
    for(int i = 0; i <= c; i++)
        arq << 0 << " ";
    arq << endl;

	//imprimir os padroes em um arquivo
    list< list<float> >::iterator it;
    for (it = padrao.begin(); it != padrao.end(); ++it){
	list<float>tl=*it;
        list<float>::iterator itt;	
	for (itt= tl.begin(); itt != tl.end(); itt++)
		arq << *itt << " ";
	arq << endl;
    }
	
	/**Modo mais simples de fazer o loop acima(WINDOWS only)
    for (list< list<float> >::iterator L : padrao) {
       for (float x : L)
          arq << x << " ";
       arq << endl;
    }*/
 
    arq.close();
    delete modelo;

    cout << "Padroes gerados com sucesso!" << endl;
    //cout << "Numero de solucoes: " << contador << endl;
    //cout << "Tempo gasto: " << time << endl;
	


    //--------------------------------------------------------CPLEX part
	int fo;
	cout << "Digite a funcao objetivo a ser usada: \n - 1: Minimizar sobra total \n - 2: Minimizar formas utilizadas no ultimo periodo\n" << endl;
	cin >> fo;
    try{
		Problema_Vigas1	Prob(inst, "padroes.txt");       // pelo cmd
		
		//Prob.heuristica();

        Prob.iniciar_lp(fo);

        cout << "\n\n\nResolvendo relaxacao linear... \n\n\n";	//resolver relaxacao linear
    	timeused(NULL);
		Prob.resolver_linear();
    	timeused(&time2);
        cout << "\n\nTempo para resolver relaxacao linear: "<< time2 <<"\n\n";
        Prob.imprimir_solucao();
    }
    catch(...){
        cerr << endl << "\n Erro na resolucao da linear" << endl;
    }

    try {
		Problema_Vigas1	Prob(inst, "padroes.txt");
		cout << "\n\n\nResolvendo Inteira... \n\n";
        Prob.iniciar_lp(fo);
        Prob.exportar_lp();                   //criar arquivo .lp

        timeused(NULL);
        Prob.revolver_ppl();                    //resolver problema
        timeused(&time2);

        cout << "\n\nTempo do gecode + resolucao do CPLEX gasto (Solucao Inteira): " << time1 + time2 << endl;
        Prob.imprimir_solucao();
    }
    catch (...) {
        cerr << endl << "\n Erro na resolucao da inteira" << endl;
    }

    //int saida = system("pause");
    return 0;
}
