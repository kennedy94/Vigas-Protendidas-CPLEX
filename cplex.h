
#include <ilcplex/ilocplex.h>
#include "gecode.h"
#include <algorithm>
#include <iomanip>
#include <vector>


class Padrao{
public:
    int tipo;
    int *tamanhos;
    float cap;
    //int E;

	bool contem(int tam){
		if(tamanhos[tam] > 0)
			return true;
		return false;
	}
};

class Solucao{

protected:
	vector< vector<Padrao> > _PERIODOS;
	int M, T;
public:
	vector< vector<Padrao> > get_Periodos(){
		return _PERIODOS;
	}

	int get_T(){
		return _PERIODOS.size();
	}

	Solucao(int _T, int _M, Padrao Pzero){
		T = _T; M = _M;
		_PERIODOS.resize(T);
		for(int t = 0; t < T; t++)
			_PERIODOS[t].resize(M);

		for(int t = 0; t < T; t++)
			for(int m = 0; m < M; m++)
				_PERIODOS[t][m] = Pzero;

	}

	void adicionar_periodo(Padrao P, int t, int m){
		_PERIODOS[t][m] = P;
	}

	void imprimir(int novo_T){
		T = novo_T;
		_PERIODOS.resize(T);
		for(int t = 0; t < T; t++){
			for(int m = 0; m < M; m++)
				cout <<  _PERIODOS[t][m].cap << " ";
			cout << endl;
		}
	}
};

class Problema_Vigas1 {
	
private:
    //dados da instancia
    int C, M, T;
    float *c_;

    Tipo_Viga *Viga;
    int P; //#padroes
    Padrao *Pattern;
    int ***X;
    //variaveis do cplex
    IloEnv env;
    IloModel model;
    IloCplex cplex;
    IloBoolVarArray x;          //1 se o padrao 'i' comeca na forma 'm' no dia 't'
	//IloBoolVarArray z;
	bool relaxacaolinear;
	int numvar, numvarP;
public:

	Padrao get_Padroes(){
		return *Pattern;
	}

    //construtor para ler os arquivos
    Problema_Vigas1(const char* filename, const char* filename2) {

        //leitura da instancia ---------------------------------------------------
        ifstream instancia(filename, ifstream::in);
        if (instancia.fail()) {
            cerr << "     Arquivo \"" << filename << "\" nao encontrado." << endl;
            exit(1);
        }
        instancia >> C >> M >> T;
        c_ = new float[M];
        Viga = new Tipo_Viga[C];
        for(int i = 0; i < C; ++i){
            instancia >> Viga[i].e >> Viga[i].k;
            Viga[i].l = new float[Viga[i].k];
            Viga[i].d = new int[Viga[i].k];
            for (int j = 0; j < Viga[i].k; j++)
                instancia >> Viga[i].l[j];
            for (int j = 0; j < Viga[i].k; j++)
                instancia >> Viga[i].d[j];
        }
        for (int i = 0; i < M; i++)
            instancia >> c_[i];
        instancia.close();
 
        //leitura dos padroes
        ifstream padroes(filename2, ifstream::in);
        if (padroes.fail()) {
            cerr << "     Arquivo \"" << filename2 << "\" nao encontrado." << endl;
            exit(1);
        }
        padroes >> P;
        Pattern = new Padrao[P];
         
        for(int i = 0; i < P; i++){
            padroes >> Pattern[i].tipo;
            	//cout << Pattern[i].tipo<< " ";
		//cout << Viga[Pattern[i].tipo].k;
            Pattern[i].tamanhos = new int[Viga[Pattern[i].tipo].k];
            //Pattern[i].E = 0;
            for(int j = 0; j < Viga[Pattern[i].tipo].k; j++){
                padroes >> Pattern[i].tamanhos[j];
                //Pattern[i].E =+ Viga[Pattern[i].tipo].e*Pattern[i].tamanhos[j];           
                //cout << Pattern[i].tamanhos[j] << " ";
                 
            }
            padroes >> Pattern[i].cap;
            //cout << Pattern[i].cap << endl;
			//if(Pattern[i].tipo == 1)
			//	system("pause");
        }	
        padroes.close(); 
    }
	bool checar_demanda(Solucao S, int c, int k){
		vector< vector<Padrao> > _PERIODOS = S.get_Periodos();
		int soma = 0;
		for(int t = 0; t < T; t++)
			for(int m = 0; m < M; m++)
				soma += _PERIODOS[t][m].tamanhos[k];
		if (soma >= Viga[c].d[k])
			return true;
		return false;
	}
    void heuristica(){
		Solucao solu(T, M, Pattern[0]);

		int t = 0;
		int m = 0;
		
		for (int c = 0; c < C; c++)
		{
			for (int k = 0; k < Viga[c].k; k++)
			{
				for (int i = 0; i < P; i++)
				{
					if(maximal(Pattern[i], c_[m]))
					{
						if(checar_demanda(solu, c, k))
							break;
						if(Pattern[i].contem(k)){
							
							for(int cura = 0; cura < Viga[c].e; cura++)
								solu.adicionar_periodo(Pattern[i], t+cura, m);

							m++;
							
						}
						else
							continue;
						
						if(m == M){
							m = 0;
							t++;
						}
						if(t == T){
							cout << "T dado melhor que calculado pela heuristica" << endl;
							return;
						}
					}
				}
			}
		}
		solu.imprimir(t);
    }
 
    void iniciar_variaveis() {
		//z = IloBoolVarArray(env, T);

        X = new int**[P];
        for (int i = 0; i < P; i++) {
            X[i] = new int*[M];
            for (int m = 0; m < M; m++) {
                X[i][m] = new int[T];
            }
        }
          
        int contador = 0;
        for (int i = 0; i < P; i++){
            for (int m = 0; m < M; m++){
                for (int t = 0; t < T; t++){
					if(maximal(Pattern[i], c_[m]) || i == 0){
						X[i][m][t] = contador++;
					}
                }
            }
        }
        numvar = P*M*T;
        numvarP = contador + 1;


		//system("pause");
        x = IloBoolVarArray(env, contador+1);//P*M*T);
        int h = 0;
        char strnum[30];

        //nome bonitinho para variavel
        for (int i = 0; i < P; i++){
            for (int m = 0; m < M; m++){
                for (int t = 0; t < T; t++){
					if(maximal(Pattern[i], c_[m]) || i == 0){
						sprintf_s(strnum, "x(%d,%d,%d)", i, m, t);
						x[h].setName(strnum);
						++h;
					}
                }
            }
        }
 
    }
    //lembrete: c_ eh um vetor de float
    void funcao_objetivo() {
        IloInt m, i, t;
  
        IloExpr costSum(env);
         
        for (m = 0; m < M; m++) 
            for (t = 0; t < T; t++)
                for (i = 1; i < P; i++)
                    if(Pattern[i].cap <= c_[m] && maximal(Pattern[i], c_[m]))
                        costSum += Viga[Pattern[i].tipo].e * (c_[m]-Pattern[i].cap) * x[X[i][m][t]];
 
        model.add(IloMinimize(env, costSum)).setName("FO");
        costSum.end();
    }


	void funcao_objetivo2(){
		IloInt m, j, i;


		IloExpr costSum(env);
		for( m = 0; m < M; m++)
			for(i = 0; i < P; i++)
				if(	maximal(Pattern[i], c_[m]) && (Pattern[i].cap <= c_[m]) || i == 0)
					costSum += x[X[i][m][T-1]];


		model.add(IloMinimize(env, costSum)).setName("FO#2");
        costSum.end();

	}
 

    void restricoes_onlyone() {
        IloInt m, t, i;
        //para cada forma m e periodo de tempo t
        for (m = 0; m < M; m++){
            for (t = 0; t < T; t++){
                IloExpr expr(env);
                expr += x[X[0][m][t]];
                for (i = 1; i < P; i++){
                    if ((Pattern[i].cap <= c_[m]) && (t + Viga[Pattern[i].tipo].e <= T) && maximal(Pattern[i], c_[m])) {
                            expr += x[X[i][m][t]];
                    }
                }
                model.add(expr <= 1).setName("Um Padrao");
                expr.end();
            }
        }
    }
    void restricoes_demanda(){
 
        IloInt c, k, m, t, i;
        //Para cada tipo de viga e tamanhos dentro do tipo de viga
        for(c = 0; c < C; c++){
            for(k = 0; k < Viga[c].k; k++){
                //sum
                IloExpr expr(env);
				//problema: alguns problemas não estão gerando padrões com a última demanda
                for(m = 0; m < M; ++m){
                    for(i = 1; i < P; ++i){
						for(t = 0; t < T - Viga[Pattern[i].tipo].e +1 ; ++t){
                            if(Pattern[i].cap <= c_[m]   &&  Pattern[i].tipo == c && maximal(Pattern[i], c_[m])){
                                expr += Pattern[i].tamanhos[k] * x[ X[i][m][t]  ];
							}
                        }
                    }
                }
                model.add(expr >= Viga[c].d[k]).setName("Demanda");
                expr.end();
            }
        }
    }
 
    void restricoes_sequenciamento(){
        int i, j, m, a, t;
        //para cada forma m, padrao t que cabe na forma m e periodo de tempo t.
		/**
        for(m = 0; m < M; ++m){
            for(i = 1; i < P; ++i){
                if(Pattern[i].cap <= c_[m] && maximal(Pattern[i], c_[m]) || i == 0){
                    for (t = 0; t < T - Viga[Pattern[i].tipo].e + 1 ; t++){                 
                        IloExpr expr(env);
                        for(j = 1; j < P; ++j)
                            if(Pattern[j].cap <= c_[m] && (maximal(Pattern[j], c_[m]) || j == 0))
                                for(a = 1; a <= Viga[Pattern[i].tipo].e - 1 ; a++)
									expr += x[X[j][m][t+a]];

                        expr += Viga[Pattern[i].tipo].e * x[    X[i][m][t]  ];
                        model.add(expr <= Viga[Pattern[i].tipo].e).setName("Sequenciamento");
                        expr.end();
                    }
                }
            }
        }
		
		*/
		for(m = 0; m < M; ++m)
			for(i = 1; i < P; ++i)
                if(Pattern[i].cap <= c_[m] && maximal(Pattern[i], c_[m]))
					for (t = 0; t < T - Viga[Pattern[i].tipo].e + 1 ; t++){
						IloExpr expr(env);
						for(a = 1; a <= Viga[Pattern[i].tipo].e - 1; a++)
							expr += x[X[0][m][t+a]];
						if((Viga[Pattern[i].tipo].e - 1) != 0)
							model.add((Viga[Pattern[i].tipo].e - 1) * x[X[i][m][t]] <=  expr).setName("problema");
						expr.end();
					}

		
		int R = 0;
		for(int c = 0; c < C; c++)
			if(Viga[c].e > R)
				R = Viga[c].e;
		
		for(t = 0; t < T; t++)
			for(m = 0; m < M; m++){
				IloExpr expr(env);

				for(int beta = 2; beta <= R; beta++)
					for(j = beta; j <= R; j++)
						for(i = 0; i < P; i++)
							if(Pattern[i].cap <= c_[m] && maximal(Pattern[i], c_[m]) && Viga[Pattern[i].tipo].e == j && t-beta+1 >= 0)
								expr += x[X[i][m][t-beta+1]];

				model.add(x[X[0][m][t]] <= expr).setName("oi");
				expr.end();
			}
    }
	
	void restricoes_continuidade(){
		IloInt m, t, i;
		
		for(m = 0; m < M; m++){
			for(t = 0; t < T-1; t++){
				IloExpr expr1(env), expr2(env);

				for(i = 0; i < P; i++)
					if(maximal(Pattern[i], c_[m]) && (Pattern[i].cap <= c_[m]) || i == 0)
						expr1 += x[X[i][m][t]];
				for(i = 0; i < P; i++)
					if(maximal(Pattern[i], c_[m]) && (Pattern[i].cap <= c_[m]) || i == 0)
						expr2 += x[X[i][m][t+1]] * Pattern[i].cap;

				model.add(c_[m] * expr1 >= expr2);
				expr1.end();
				expr2.end();
			}
		}
		
	}
 
    void exportar_lp() {
        cplex.exportModel("problema_vigas.lp");
    }

	void resolver_linear(){
		relaxacaolinear = true;
        IloModel relax(env);
        relax.add(model);
        relax.add(IloConversion(env, x, ILOFLOAT)); 
        cplex = IloCplex(relax);
        if (!cplex.solve()) {
            env.error() << "Otimizacao do LP mal-sucedida." << endl;
            throw(-1);
        }
    }
    void revolver_ppl() {
	//cplex.setParam(IloCplex::PreInd, 0); Desligar presolve(NAO FACA ISSO DE NOVO!)

        cout << "Variaveis sem preprocessamento: " << numvar << endl;
    	cout << "Variaveis com preprocessamento: " << numvarP << endl << endl;
    	cplex.setParam(IloCplex::TiLim, 600);
    	//cplex.setParam(IloCplex::Param::MIP::Cuts::Cliques, -1);
            if (!cplex.solve()) {
                env.error() << "Otimizacao do LP mal-sucedida." << endl;
                throw(-1);
            }
    }
    void imprimir_solucao() {
  
        cplex.out() << "Status da solucao = " << cplex.getStatus() << endl;
        cplex.out() << "Valor Otimo  = " << cplex.getObjValue() << endl;
        cplex.out() << "#Iteracoes = " << cplex.getNiterations() << endl;
        cplex.out() << "#Nos de BB  = " << cplex.getNnodes() << endl;
        cout << "\n     Funcao Objetivo: " << cplex.getObjValue() << endl;
		
		bool imprimir = 1;
		if(imprimir){
			int contador = 1;
			for (int t = 0; t < T; t++) {
				for (int m = 0; m < M; m++) {
					for (int i = 0; i < P; i++) {
                
							if(maximal(Pattern[i], c_[m]) || i == 0){
								if(cplex.isExtracted(x[X[i][m][t]]) && cplex.getValue(x[X[i][m][t]]) > 0.00001){
									cout << "     x(" << i << "," << m << "," << t << ") = " << cplex.getValue(x[X[i][m][t]]) <<"		";
									if(contador%2 == 0)
										cout << endl;
									contador++;
								}
							}
						}
				}
			}   
			cout << "\n\n";
			//se for relaxação linear não imprimir Gantt
			if(relaxacaolinear)
				return;

			//cout << "Variaveis sem preprocessamento: " << numvar << endl;
			//cout << "Variaveis com preprocessamento: " << numvarP << endl << endl;

			// usando a biblioteca "#include <iomanip>"
			cout << "  \"Gantt\" " << endl << endl;
			const char separator    = ' ';
			const char separator2    = '_';
			const int nameWidth     = 6;
			//const int numWidth      = 8;

			for (int t = 0; t <= T; t++)
				cout << internal << setw(nameWidth) << setfill(separator2) << separator2 ;
			cout << endl;
			cout << internal << setw(nameWidth) << setfill(separator) << separator;
			for (int t = 0; t < T; t++)
				cout << internal << setw(nameWidth) << setfill(separator) << t ;
		

			//imprimir sobra = 1, padrao = 0, tipo do padrao = 2?
			int gantt = 2	;
			//--------------
			cout << endl;


			for (int m = 0; m < M; m++) {
				cout << internal << setw(nameWidth) << setfill(separator) << m;
				for (int t = 0; t < T; t++) {
					bool nenhum = 1;
					for (int i = 0; i < P; i++) {
						if((maximal(Pattern[i], c_[m]) || i == 0) && cplex.isExtracted(x[X[i][m][t]]) && cplex.getValue(x[X[i][m][t]]) > 0.00001){
							nenhum = 0;
							if(i == 0){
								cout << internal << setw(nameWidth) << setfill(separator) << "*";
								break;
							}
							else{

								switch (gantt)
								{
								case 0:
									cout << internal << setw(nameWidth) << setfill(separator) << i; //imprimindo numero do padrao
									break;
								case 1:
									 //imprimindo sobra
									if(i == 0)
										cout << internal << setw(nameWidth) << setfill(separator) << "*";
									else
										cout << internal << setw(nameWidth) << setfill(separator) << c_[m] - Pattern[i].cap;
									break;
								case 2:
									if(i == 0)
										cout << internal << setw(nameWidth) << setfill(separator) << "*";
									else
										cout << setprecision(2) << internal << setw(nameWidth) << setfill(separator) << Pattern[i].cap/c_[m] ;
									break;
								}
							}		
						}	
					}
					if(nenhum)
								cout << internal << setw(nameWidth) << setfill(separator) << " ";
				}
				cout << endl;
			}
		
			for (int t = 0; t <= T; t++)
				cout << internal << setw(nameWidth) << setfill(separator2) << separator2 ;
			cout << endl;

		}

		if(verificacao())
			cout << "\n\nSolucao valida para instancia!\n\n";
		else
			cout << "\n\nSolucao nao válida para instancia!\n\n";
    }


	//----------------------------------------------------------------------
	bool maximal(const Padrao P, float C_FORMA){
        //return true;
		float menor = 10000;
		for(int aux = 0; aux < Viga[P.tipo].k; aux++)
			if (Viga[P.tipo].l[aux] < menor) menor = Viga[P.tipo].l[aux];
        //cout << "------" <<C_FORMA << " - " << P.cap << " < ? " << menor << endl;
 		return (C_FORMA - P.cap) < menor; 

	}/**
    void verificar_maximais(){

        for(int i = 1; i < P; i++)
            for(int m = 0; m < M; m++)
                cout << "Padrao " << i << " e maximal na forma " << m << "?" << maximal(Pattern[i], c_[m]) << endl;

    }*/
	//pseudo-funcao de verificacao
	bool verificacao(){
	//----------------------------------------------------------------------
		//obedece sequenciamento?
		int h = 0;
		for(int m = 0; m < M; m++){
			for(int t = 0; t < T; t++){
				for (int i = 1; i < P; i++) {
					if(maximal(Pattern[i], c_[m])){
						if(cplex.isExtracted(x[X[i][m][t]]) && cplex.getValue(x[X[i][m][t]]) > 0.00001){

							if( Viga[Pattern[i].tipo].e  > 1){
								for(int a = 1; a < Viga[Pattern[i].tipo].e - 1 ; t++)
									if( cplex.getValue(x[X[0][m][t+a]]) != 1)
										return false;
							}
						}
					}
				}
			}
		}

		//apenas um padrão??
		for(int m = 0; m < M; m++){
			for(int t = 0; t < T; t++){
				for (int i = 0; i < P; i++) {
					int soma = 0;
					if(maximal(Pattern[i], c_[m]) || i == 0)
						if(cplex.isExtracted(x[X[i][m][t]]) && cplex.getValue(x[X[i][m][t]]) > 0.00001)
							soma++;

					if (soma > 1)
						return false;
				}
			}
		}

		return true;
	}
    void iniciar_lp(int fo) {
        try {
			relaxacaolinear = false;

            model = IloModel(env);

            iniciar_variaveis();

			switch (fo)
			{
			case 1:
				funcao_objetivo();
				break;
			case 2:
				funcao_objetivo2();
				break;
			}


            restricoes_onlyone();

            restricoes_demanda();

            restricoes_sequenciamento();

			//restricoes_continuidade();
			
			//imprimir quantas casas decimais?
			//int decimal = 2;
			//cout << fixed;			//imprimindo apenas 2 casas decimais
			//cout << setprecision(decimal);
	
            cplex = IloCplex(model);
        }
        catch (IloException& e) {
            cerr << "Erro: " << e.getMessage() << endl;
            cout << "\nErro ilocplex" << endl;
            return;
        }
        catch (...) {
            cerr << "Outra excecao" << endl;
            cout << "\noutra excecao" << endl;
        }
    }
  
    ~Problema_Vigas1()
    {
        env.end();/**
        //destroi tudo
        delete []l;
        delete []c;
        delete []e;
 
        for(int i = 0; i < m; i++){
            for(int j= 0; j < n; j++){
                delete []X[i][j];
            }
        }*/
    }
};
