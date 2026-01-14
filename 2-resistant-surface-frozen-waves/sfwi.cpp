#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <cmath>
#include <iomanip>      // Formatar cout
#include <ctime>        // Contar tempo decorrido
#include <cstdlib>      // Converter string em double

// VERSOES (funcao info() contem mais informacoes sobre o codigo)
// v1_6
// - otimizacao de calculos
// v1_5
// - otimizacao para P<pxH e P>pxH
// - calculo dos coeficientes Aqp's
// - funcao binarizada f(x,z) como entrada
// V1_4
// - foi retirado o indice de refracao, as sfws serao consideradas no vacuo;
// - os parametros da sfw sao passados como argumentos de linha de comando;
// - exclusao da funcao progresso;
// V1_3
// - o programa recebe como parâmetro (argumento) da linha de comando o
//caminho do arquivo MTX com os coeficientes Aqp's;
// - o arquivo de saida contendo o calculo de intensidade sera exportado como
// '[diretorio_do_arquivo_de_entrada\nome_do_arquivo_de_entrada.mtx]_SFWI.mtx';
// V1_2
// - exclusao de variaveis desnecessarias;
// - otimizacao do calculo para um melhor desempenho;

#define NOME "SFWI"
#define DATA  "16 ABR 2020"
#define VERSAO "1_6"
#define AUTOR "@JODESARRO"
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
#endif

using namespace std;

void info();

int main(int argc, char *argv[])
{

    // === DECLARAÇOES ====================================================
    // Todas as variaveis utilizadas estao declaradas aqui no comeco com
    //a intencao de prezar pelo processamento ao custo de memoria
    // Gerais
    const double M_2PI = 2.0*M_PI;// = 2*M_PI
    long int tempo = time(NULL);  // inteiro para calcular o tempo decorrido
    string file_in;               // Caminho do arquivo de entrada (funcao F em MTX)
    string file_out;              // Caminho do arquivo de saida (intensidades)
    // Frozen Waves
    double * F;                 // Array com o padrao discretizado F(x,z) (imagem binarizada)
    complex<double> * Aqp;      // Coeficientes complexos Aqp
    string tmp_s;               // Variavel temporaria tipo string
    complex<double> tmp_cd;     // Variavel temporaria tipo complex double
    double l0;                  // Lambda0, comprimento de onda
    double L;                   // Parâmetro longitudinal das FWs referente a largura da imagem gerada (em z)
    double H;                   // Parâmetro radial das FWs referente a altura da imagem gerada (em x)
    int N;                      // Parametro de feixes de Bessel para um total de 2N+1
    int N2M1;                   // Total de feixes de Bessel = N*2+1
    int P;                      // Quantidade de LFWs
    double HdP;                 // = H/P, altura que cada lfw ocuparia em x
    double Q;                   // Parâmetro Q das FWs
    double kk;                  // = k², vetor de onda ao quadrado
    double * b;                 // beta, vetor de onda longitudinal
    double * h;                 // krho, vetor de onda radial
    double * x0;                // Deslocamento x0 das LFWs
    double * x0x0;              // = x0²
    int pxL;                    // Total de pixels em z da imagem de entrada
    int pxLm1;                  // = pxL-1
    int pxH;                    // Total de pixels em x da imagem de entrada
    int pxHdP;                  // = int ceil(pxH/P)
    int pxHL;                   // = pxL*pxH, area em pixels da imagem de entrada
    double * sfwi;              // Array para armazenar o calculo da intensidade da sfw
    int pxHipdP;                // Variavel que relaciona os pixels em x da imagem de entrada com a quantidade p de LFW
    // Variaveis para calculos
    int DX;                     // Quantidade de divisoes da variavel x
    int DXm1;                   // = DX-1
    double dx;                  // Incremento da variavel x
    double var_x;               // Variavel x
    double var_xx;              // = x*x
    double var_2x;              // = 2*x
    int DZ;                     // Quantidade de divisoes da variavel z
    double dz;                  // Incremento da variavel z
    double var_z;               // Variavel z
    // ====================================================================


    // Mostrar informacoes sobre o programa
    info();


    // === ARGUMENTOS/PARÂMETROS/ENTRADA ==================================
    if ( argc != 8 )
    {
        puts("Erro ao processar algum parametro. Verifique se todos foram passados como argumento.");
        puts("O programa foi encerrado.");
        exit(1);
    }

    l0 = atof(argv[1]);         // Lambda0, comprimento de onda 
    Q  = atof(argv[2]);         // Parâmetro Q das frozen waves 
    N  = atoi(argv[3]);         // Parâmetro N dos feixes de Bessel
    L  = atof(argv[4]);         // Parâmetro longitudinal das FWs
    P  = atoi(argv[5]);         // Quantidade de LFWs 
    file_in = argv[6];          // Nome do arquivo que sera importado
    file_out   = argv[7];       // Nome do arquivo que sera exportado
    
    N2M1 = N*2+1; // = N*2+1 : Total de feixes de Bessel
    // ====================================================================




    // === FUNÇAO F(x,z) e PARÂMETROS DA SFW ==============================
    // Abrir arquivo .mtx que contem a funcao discretizada F 
    fstream file_mtx(file_in);
    if ( !file_mtx.is_open() )
    {
        cout << "Nao foi possivel abrir o arquivo '" << file_in << "'." << endl;
        puts("O programa foi encerrado.");
        exit(1);
    }

    // Leitura dos dados sobre a quantidade de linhas e de colunas da matriz,
    //que correspondem a quantidade de pixels na longitudinal (pxL) e na vertical (pxH)
    file_mtx >> pxH >> pxL;
    // Se a leitura nao corresponder a um inteiro, como por exemplo nos
    //cabecalho ou comentarios do arquivo MTX, a leitura passara para a proxima linha.
    while ( !file_mtx.good() )
    {
        file_mtx.clear();
        file_mtx.ignore(INT_MAX, '\n');
        file_mtx >> pxH >> pxL;
    }

    pxHL = pxH*pxL;     // Area da funcao de entrada em pixels
    pxLm1 = pxL - 1;    // = pxL-1
    H = L*pxH/pxL;      // Altura da imagem de saida proporcional a de entrada
    HdP = H/P;          // = H/P

    // Alocacao de memoria para o array F
    F = new double [ pxHL ];
    // Preenchendo o array F
    for ( int i=0; i<pxHL; i++ )
    {
        file_mtx >> F[i];
    }
    // Fechar o arquivo .mtx da funcao binarizada F
    file_mtx.close();
    // ====================================================================




    // ==== ARQUIVO DE SAIDA =============================================
    // Criando o arquivo de saida, pois, se nao for possivel, nem realizara o calculo
    file_mtx.open(file_out, fstream::out);
    if ( !file_mtx.is_open() )
    {
        cout << "Nao foi possivel exportar para o arquivo '" << file_out << "'." << endl;
        puts("O programa foi encerrado.");
        exit(1);
    }
    // ===================================================================


    puts("Calculando...");


    // === VETORES DE ONDA ================================================
    // Vetor de onda no meio ao quadrado
    kk = (M_2PI/l0)*(M_2PI/l0);

    // Alocacao de memoria para os vetores de ondas longitudinais e radiais 
    b = new double [ N2M1 ];
    h = new double [ N2M1 ];

    // Calculo dos vetores de onda longitudinais e radiais 
    // No bloco abaixo, iq=q+N, q=iq-N, [iq]->q 
    for ( int iq=0; iq<N2M1; iq++ )
    {
        b[iq] = Q + M_2PI*(iq-N)/L ;
        h[iq] = sqrt( kk - b[iq]*b[iq] );
    }
    // ===================================================================




    // === CALCULO DOS COEFICIENTES Aqps ==================================
    if ( P <= pxH ) // A[q][p]
    {
        // Alocacao de memoria para o array A[q][p] 
        Aqp = new complex<double> [ P*N2M1 ];
        // Iteracoes em p 
        for (int ip=0; ip<P; ip++)
        {
            // Variavel que relaciona os pixels em x da imagem de entrada com a quantidade p de LFW
            pxHipdP = int ( ceil(pxH*ip/P) );
            // Iteracoes em q
            for (int iq=0; iq<N2M1; iq++)
            {
                // O proximo for contem iteracoes em j que corresponde aos pixels em z.
                // Esta iteracao e uma aproximacao da integral
                //dos coeficientes Aqp atraves do metodo trapezoidal
                //para espacamentos uniformes em z. Cada j corresponde a
                //um ponto em z (ou um pixel em zj da imagem de entrada)
                tmp_cd = 0.5*( F[ pxHipdP + pxH*pxLm1 ] + F[ pxHipdP + 0 ] );
                for (int j=1; j<pxLm1; j++)
                {
                    tmp_cd +=  F[ pxHipdP + pxH*j ]*exp( complex<double>(0., -M_2PI*(iq-N)*j/pxLm1) ) ;                    
                }
                Aqp[ iq + N2M1*ip ] = tmp_cd/( double (pxLm1) );
            }
        }
    }
    else // Otimizacao para pxH < P, Aqp -> A[q][pxh]
    {
       // Alocacao de memoria para o array A[q][pxh]
        Aqp = new complex<double> [ pxH*N2M1 ];
        // Iteracoes em pxh
        for (int ipxh=0; ipxh<pxH; ipxh++)
        {
            // Iteracoes em q
            for (int iq=0; iq<N2M1; iq++)
            {
                // O proximo for comtem Iteracoes em j que corresponde aos pixels em z.
                // Esta iteracao e uma aproximacao da integral
                //dos coeficientes Aqp atraves do metodo trapezoidal
                //para espacamentos uniformes em z. Cada j corresponde a
                //um ponto em z (ou um pixel em zj da imagem de entrada)
                tmp_cd = 0.5*( F[ ipxh + pxH*pxLm1 ] + F[ ipxh + 0 ] );
                for (int j=1; j<pxLm1; j++)
                {
                    tmp_cd +=  F[ ipxh + pxH*j ]*exp( complex<double>(0., -M_2PI*(iq-N)*j/pxLm1) ) ;                    
                }
                Aqp[ iq + N2M1*ipxh ] = tmp_cd/( double (pxLm1) );
            }
        }
    }
    // ====================================================================




    // === DESLOCAMENTO DE CADA LFW COM RELAÇAO A ORIGEM =================
    // Alocacao de memoria para os valores dos deslocamentos 
    x0 = new double [ P ];
    x0x0 = new double [ P ];

    // No bloco abaixo, ip=p-1, p=ip+1, [ip]->p 
    for (int ip=0; ip<P; ip++)
    {
        // Distribuicao linear e homogênea em x de acordo com P e H.
        x0[ip] = (ip+1.0)*HdP;
        x0x0[ip] = x0[ip]*x0[ip];
    }
    // ===================================================================




    // === DEFINICOES VARIAVEIS DOS EIXOS X E Z ==========================
    DX = 2*P;       // quantidade de divisoes de x 
    DXm1 = DX-1;    // = DX - 1
    dx = H/DX;      // incremento em x 
    DZ = pxL;       // quantidade de divisoes de z 
    dz = L/DZ;      // incremento em z
    // ===================================================================




    // === CALCULO DA SURFACE FROZEN WAVE ================================
    // Alocacao de memoria para a sfw 
    sfwi = new double [ DX*DZ ];

    // Iteracoes na variavel longitudinal z 
    for (int iz=0; iz<DZ; iz++)
    {
        var_z = (iz+1.0)*dz;

        // Iteracoes na variavel x perpendicular a z 
        for ( int ix=0; ix<DX; ix++ )
        {
            var_x = (ix+1.0)*dx;
            var_xx = var_x*var_x;
            var_2x = 2.0*var_x;
            tmp_cd = complex<double> (0.0, 0.0);

            // Iteracoes em p 
            for (int ip=0; ip<P; ip++)
            {
                if ( P <= pxH ) // Nota: Aqp -> A[q][p]
                {
                    // Iteracoes em q 
                    for (int iq=0; iq<N2M1; iq++)
                    {
                        tmp_cd += Aqp[ iq + N2M1*ip ] * cyl_bessel_j( 0, h[iq]*sqrt(var_xx+x0x0[ip]-var_2x*x0[ip]) ) * exp( complex<double> (0., b[iq]*var_z) );
                    }
                }
                else // Devido a otimizacao para pxH < P, Aqp -> A[q][pxh]
                {
                    pxHipdP = int ( ceil(pxH*ip/P) );
                    // Iteracoes em q 
                    for (int iq=0; iq<N2M1; iq++)
                    {
                        tmp_cd += Aqp[ iq + N2M1*pxHipdP ] * cyl_bessel_j( 0, h[iq]*sqrt(var_xx+x0x0[ip]-var_2x*x0[ip]) ) * exp( complex<double> (0., b[iq]*var_z) );
                    }
                }
            }
            // Modulo ao quadrado da funcao de onda
            sfwi[ix + DX*iz] = abs(tmp_cd)*abs(tmp_cd);
        }
    }
    // ===================================================================




    // === EXPORTAÇAO DA INTENSIDADE, |sfw|², DA SFW =====================
    // O arquivo de saida ja fora criado antes do calculo
    // Escrita dos cabecalho e comentario no arquivo MTX
    //que corresponde às suas duas primeiras linhas      
    file_mtx << "%%MatrixMarket matrix array real general\n";                  // Cabecalho do mtx 
    file_mtx << "%Criado com o software " << NOME << " " << VERSAO << endl;    // Primeiro comentario do mtx  
    
    // Escrita de dados da matriz na terceira linha do arquivo 
    file_mtx << DX << " " << DZ << endl;

    // Iteracoes na variavel longitudinal z 
    for (int iz=0; iz<DZ; iz++)
    {
        var_z = (iz+1.0)*dz;

        // Iteracoes na variavel x perpendicular a z 
        for ( int ix=0; ix<DX; ix++ )
        {
            var_x = (ix+1.0)*dx;
            file_mtx << setprecision(16) << scientific << "   " << sfwi[ix + DX*iz] << endl;
        }
    }

    // Fechar o arquivo .mtx dos coeficientes Aqp 
    file_mtx.close();
    // ===================================================================




    // === LIBERAÇAO DE MEMÓRIAS ALOCADAS ================================
    delete[] Aqp;
    delete[] x0;
    delete[] x0x0;
    delete[] b;
    delete[] h;
    delete[] sfwi;
    // ===================================================================




    // === FINALIZANDO ===================================================
    tempo = time(NULL) - tempo;
    cout << "O processo levou " << tempo << " segundos." << endl;
    // ===================================================================


}

void info()
{
    cout << "| " << NOME << " " << VERSAO << " | " << DATA << " | " << AUTOR << " |" << endl;
    puts("");
    puts("DESCRICAO");
    puts("Calculo de intensidade de Surface Frozen Wave (SFW) ordinaria no plano 'xz', no vacuo, com altura e largura proporcionais a funcao F(x,z) e LFWs igualmente espacadas no eixo 'x' e de mesmo parametros 'L' e 'N' ,");
    puts("");
    puts("ENTRADA");
    puts("Os seguintes parametros devem ser passados como argumentos ao programa, separados por espacos, respectivamente:");
    puts(" - Comprimento de onda no vacuo em metros;");
    puts(" - Parametro 'Q' das Frozen Waves;");
    puts(" - Parametro 'N' dos feixes de Bessel para um total de 2N+1;");
    puts(" - Parametro longitudinal 'L' (largura) da SFW em metros;");
    puts(" - Quantidade 'P' de LFWs;");
    puts(" - Caminho do arquivo no formato array Matrix Market (MTX) contendo a funcao discretizada F(x,z);");
    puts(" - Caminho do arquivo que sera criado pelo programa para armazenar a intensidade calculada da SFW em formato MTX;");
    puts("Exemplo de execucao: 'sfwi.exe' .000000632 9939761. 5 .06 328 'aqp.mtx' 'int.mtx'"); 
    puts("");
    puts("SAIDA");
    puts("Arquivo no formato MTX contendo as intensidades (modulo ao quadrado) da SFW.");
    puts("");
    puts("PROGRESSO");
}