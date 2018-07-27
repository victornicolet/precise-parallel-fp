/* File for sequential viterbi implementations.
 * Author: Raphaël Dang-Nhu
 * Date: July 26th */

#include <iostream>
#include <emmintrin.h>
#include <mpreal.h>
#include <mpfr.h>
#include <gmp.h>

#include "interval_arithmetic.hpp"
#include "debug.hpp"

using mpfr::mpreal;


void viterbi_double(long n, double** t){

	double** p = new double*[n];
    for(int j = 0; j != n; j++){
        p[j] = new double[n];
    }

	for(int i = 0; i != n; i++){
		p[0][i] = 1.;
	}

	for(int i = 1; i != n; i++){
		for(int j = 0; j != n; j++){
			p[i][j] = 0.;

			for(int k = 0; k != n; k++){
				
				double aux = p[i-1][k]*t[k][j];

				if(aux > p[i][j]){
					p[i][j] = aux;
				}

			}
		}
	}
	
	double result = 0.;
	for(int j = 0; j != n; j++){
		if(p[n-1][j] > result){
			result = p[n-1][j];
		}
	}

    if(PRINT){
        cout << endl << "Viterbi" << endl;
        cout << "Result: " << result << endl;
    }

}

void viterbi_mpfr( long n, double** t){

	mpreal** p = new mpreal*[n];
    for(int j = 0; j != n; j++){
        p[j] = new mpreal[n];
    }

	for(int i = 0; i != n; i++){
		p[0][i] = mpreal(1.,1000);
	}

	for(int i = 1; i != n; i++){
		for(int j = 0; j != n; j++){
			p[i][j] = mpreal(0.,1000);

			for(int k = 0; k != n; k++){
				
				mpreal aux = p[i-1][k]*t[k][j];

				if(aux > p[i][j]){
					p[i][j] = aux;
				}

			}
		}
	}
	
	mpreal result = 0.;
	for(int j = 0; j != n; j++){
		if(p[n-1][j] > result){
			result = p[n-1][j];
		}
	}

    if(PRINT){
        cout << endl << "Viterbi Mpfr" << endl;
        cout << "Result: " << result << endl;
    }

}

void viterbi_interval( long n, double** t,boolean*& m){
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    long size = 2*n+1+((1 + 2*n) * n) * (n - 1);
    m = new boolean[size];
    long it = 0;
    boolean b;
    
    //Interval Approximation
	__m128d** p = new __m128d*[n];
    for(int j = 0; j != n; j++){
        p[j] = new __m128d[n];
    }
    __m128d aux;
    __m128d result;
    
    for(int i = 0; i < n;i++){
        p[0][i]= in2_create(1.);
        it++;  
    }
    
    for(int i___0 = 1; i___0 < n;i___0++){
        
        for(int j = 0; j < n;j++){
            p[i___0][j] = in2_create(0.);
            it++; 
            for(int k = 0; k < n;k++){
                aux = in2_mul(p[i___0 - 1][k],t[k][j]);
                it++; 
                b = in2_gt(aux,p[i___0][j]);
                m[it] = b;
                it++;
                
                switch(b){
                    case True:
                        p[i___0][j] = aux; 
                    break;
                    case False:
                        
                    break;
                    case Undefined:
                        p[i___0][j] = in2_merge(aux,p[i___0][j]); 
                        
                    break;
                }
            }
        }
    }
    result = in2_create(0.);
    it++; 
    for(int j___0 = 0; j___0 < n;j___0++){
        b = in2_gt(p[n - 1][j___0],result);
        m[it] = b;
        it++;
        
        switch(b){
            case True:
                result = p[n-1][j___0];
            break;
            case False:
                
            break;
            case Undefined:
                result = in2_merge(p[n-1][j___0],result);
            break;
        }
         
    }

    if(PRINT){
        cout << endl << "Interval viterbi"  << endl;
        cout << "Result:";
        print(result);
        cout << endl;
        /*for(int i = 0;  i != size; i++){
            Dec d = da[i]; 
            if(d == True) cout << "X";
            else if(d == False) cout << "_";
            else cout << "|";
        }*/
        cout << endl;
    }
    _MM_SET_ROUNDING_MODE(0);
}

void viterbi_reverse( long n, double** t,boolean*m){

    long it = 2*n+1+((1 + 2*n) * n) * (n - 1);

	bool** p_rev = new bool*[n];
    for(int j = 0; j != n; j++){
        p_rev[j] = new bool[n]();
    }

    boolean b;
    bool aux_rev = false;
    bool result_rev = true;;

    for(int j___0_rev = n-1; j___0_rev >= 0;j___0_rev--){
        it--;
        b = m[it]; 
        switch(b){
            case True:
                if(! result_rev){
                    m[it] = Useless;  
                }
                if(result_rev){
                    result_rev = 0;  
                    p_rev[n - 1][j___0_rev] = 1.;
                }
            break;
            case False:
                    m[it] = Useless;  
                
            break;
            case Undefined:
                if(! result_rev){
                    m[it] = Useless;  
                }
                // Variable duplication
                
                
                
                // Left branch
                if(result_rev){
                    result_rev = 0;  
                    p_rev[n - 1][j___0_rev] = 1.;
                }
                
                //Right branch
                
                
                //Merging
                
            break;
        } 
    }
    it--; 
    if(result_rev){
        result_rev = 0;  
        m[it] = True;
    }else{
        m[it]= False;
    }
    for(int i___0_rev = n-1; i___0_rev >= 1;i___0_rev--){
        for(int j_rev = n-1; j_rev >= 0;j_rev--){
            for(int k_rev = n-1; k_rev >= 0;k_rev--){
                it = (it + -1);
                b = *(m + it); 
                switch(b){
                    case True:
                        if(! p_rev[i___0_rev][j_rev]){
                            *(m + it) = Useless;  
                        }
                        if(*(*(p_rev + i___0_rev) + j_rev)){
                            *(*(p_rev + i___0_rev) + j_rev) = 0;
                            aux_rev = 1;  
                        }
                    break;
                    case False:
                        if(! 0){
                            *(m + it) = Useless;  
                        }
                        
                    break;
                    case Undefined:
                        if(! p_rev[i___0_rev][j_rev]){
                            *(m + it) = Useless;  
                        }

                        if(*(*(p_rev + i___0_rev) + j_rev)){
                            aux_rev = 1;  
                        }
                        
                    break;
                }
                it = (it + -1); 
                if(aux_rev){
                    aux_rev = 0;
                    *(*(p_rev + (i___0_rev - 1)) + k_rev) = 1;  
                    m[it] = True;
                }else{
                    m[it]= False;
                }
            }
            it = (it + -1); 
            if(*(*(p_rev + i___0_rev) + j_rev)){
                *(*(p_rev + i___0_rev) + j_rev) = 0;  
                m[it]= True;
            }else{
                m[it]= False;
            }
        } 
    }
    for(int i_rev = n-1; i_rev >= 0;i_rev--){
        it = (it + -1); 
        if(*(*(p_rev + 0) + i_rev)){
            *(*(p_rev + 0) + i_rev) = 0;  
            m[it] = True;
        }else{
            m[it]= False;
        }
    }
    if(PRINT){
        cout << endl << "Reverse Viterbi Finished" << endl;
    }
}

void viterbi_lazy( long n, double** t,boolean*m){

    long it = 0;

	mpreal** p_ex = new mpreal*[n];
    for(int j = 0; j != n; j++){
        p_ex[j] = new mpreal[n];
    }

    boolean b;
    mpreal aux_ex (0.,1000);
    mpreal result_ex (0.,1000);
    for(int i_ex = 0; i_ex < n;i_ex++){
        b = *(m + it);
        it = (it + 1); 
        p_ex[0][i_ex] = mpreal(1.,1000);  
    }
    for(int i___0_ex = 1; i___0_ex < n;i___0_ex++){
        for(int j_ex = 0; j_ex < n;j_ex++){
            b = *(m + it);
            it = (it + 1); 
            if(b == True){
                p_ex[i___0_ex][j_ex] = mpreal(0.,1000);  
            }
            for(int k_ex = 0; k_ex < n;k_ex++){
                b = *(m + it);
                it = (it + 1); 
                if(b == True){
                    //cout << i___0_ex << "," << j_ex << endl;
                    aux_ex = p_ex[i___0_ex-1][k_ex] * t[k_ex][j_ex];  
                }
                b = *(m + it);
                it = (it + 1); 
                switch(b){
                  case True:
                        p_ex[i___0_ex][j_ex] = aux_ex; 
                    break;
                    case False:
                        
                    break;
                    case Undefined:
                        if(aux_ex > p_ex[i___0_ex][j_ex]){
                             p_ex[i___0_ex][j_ex] = aux_ex;  
                      }
                    break;
                } 
            } 
        } 
    }
    b = m[it];
    it++; 
    result_ex = 0.000000;  
    for(int j___0_ex = 0; j___0_ex < n;j___0_ex++){
        b = *(m + it);
        switch(b){
            case True:
                result_ex = p_ex[n-1][j___0_ex]; 
            break;
            case False:
                
            break;
            case Undefined:
                if(p_ex[n-1][j___0_ex] > result_ex){
                    result_ex = p_ex[n-1][j___0_ex]; 
              }else{
                     
                }
            break;
        } 
        it = (it + 1); 
    }
    if(PRINT){
        cout << endl << "Lazy Viterbi" << endl;
        cout << "Result: " << result_ex << endl;
    }
}

