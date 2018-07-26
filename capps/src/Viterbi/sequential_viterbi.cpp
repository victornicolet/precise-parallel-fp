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


void viterbi_double(double*aa, long n, double** t){

	double** p = new double*[n];
    for(int j = 0; j != n; j++){
        p[j] = new double[n];
    }

	for(int i = 0; i != n; i++){
		p[i][0] = 1.;
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
		if(p[n-1][j] < result){
			result = p[n-1][j];
		}
	}

    if(PRINT){
        cout << endl << "Viterbi" << endl;
        cout << "Result: " << *result << endl;
    }

}

void viterbi_mpfr(double*aa, long n, double** t){

	mpreal** p = new mpreal*[n];
    for(int j = 0; j != n; j++){
        p[j] = new mpreal[n];
    }

	for(int i = 0; i != n; i++){
		p[i][0] = mpreal(1.,1000);
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
		if(p[n-1][j] < result){
			result = p[n-1][j];
		}
	}

    if(PRINT){
        cout << endl << "Viterbi Mpfr" << endl;
        cout << "Result: " << *result << endl;
    }

}

void viterbi_interval(double*aa, long n, double** t,boolean**da){
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    long size = 2*n+1+((1 + 2*n) * n) * (n - 1);
    da = new boolean[size];
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
        p[i][0]= in2_create(1.);
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
        b = in2_lt(p[n - 1][j___0],result);
        m[it] = b;
        it++;
        
        switch(b){
            case True:
                result = p[n-1][j__0];
            break;
            case False:
                
            break;
            case Undefined:
                result = in2_merge(p[n-1][j__0]);
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

void viterbi_reverse(double*aa, long n, double** t,boolean**m){

	bool** p_rev = new bool*[n];
    for(int j = 0; j != n; j++){
        p_rev[j] = new bool[n];
    }

    bool aux_rev;
    bool result_rev;

    for(int j___0_rev = n-1; j___0_rev >= 0;j___0_rev--){
        it--;
        b = m[it]; 
        switch(b){
            case True:
                if(! result_rev){
                    m[it] = Useless;  
                }
                if(p_rev[n - 1][j___0_rev]){
                    p_rev[n - 1][j___0_rev] = 0.;
                    result_rev = 1;  
                }
            break;
            case False:
                if(! 0){
                    m[it] = Useless;  
                }
                
            break;
            case Undefined:
                if(! result_rev){
                    m[it] = Useless;  
                }
                // Variable duplication
                
                
                
                // Left branch
                if(p_rev[n - 1][j___0_rev]){
                    p_rev[n - 1][j___0_rev] = 0;
                    result_rev = 1;  
                }
                
                //Right branch
                
                
                //Merging
                
            break;
        } 
    }
    it--; 
    if(result_rev){
        result_rev = 0;  
    }
    for(int i___0_rev = n-1; i___0_rev >= 1;i___0_rev--){
        for(int j_rev = n-1; j_rev >= 0;j_rev--){
            for(int k_rev = n-1; k_rev >= 0;k_rev--){
                it = (it + -1);
                b = *(m + it); 
                switch(b){
                    case True:
                        if(! aux_rev){
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
                        if(! aux_rev){
                            *(m + it) = Useless;  
                        }
                        // Variable duplication
                        
                        
                        
                        // Left branch
                        if(*(*(p_rev + i___0_rev) + j_rev)){
                            *(*(p_rev + i___0_rev) + j_rev) = 0;
                            aux_rev = 1;  
                        }
                        
                        //Right branch
                        
                        
                        //Merging
                        
                    break;
                }
                it = (it + -1); 
                if(aux_rev){
                    aux_rev = 0;
                    *(*(t + k_rev) + j_rev) = 1;
                    *(*(p_rev + (i___0_rev - 1)) + k_rev) = 1;  
                } 
            }
            it = (it + -1); 
            if(*(*(p_rev + i___0_rev) + j_rev)){
                *(*(p_rev + i___0_rev) + j_rev) = 0;  
            } 
        } 
    }
    for(int i_rev = n-1; i_rev >= 0;i_rev--){
        it = (it + -1); 
        if(*(*(p_rev + i_rev) + 0)){
            *(*(p_rev + i_rev) + 0) = 0;  
        } 
    }
}

void viterbi_reverse(double*aa, long n, double** t,boolean**m){
    mpreal** p_ex;
    mpreal aux_ex (0.,1000);
    mpreal result_ex (0.,1000);
    for(int i_ex = 0; i_ex < n;i_ex++){
        b = *(m + it);
        it = (it + 1); 
        if(b == True){
            *(*(p_ex + i_ex) + 0) = 1.000000;  
        } 
    }
    for(int i___0_ex = 1; i___0_ex < n;i___0_ex++){
        for(int j_ex = 0; j_ex < n;j_ex++){
            b = *(m + it);
            it = (it + 1); 
            if(b == True){
                *(*(p_ex + i___0_ex) + j_ex) = 0.000000;  
            }
            for(int k_ex = 0; k_ex < n;k_ex++){
                b = *(m + it);
                it = (it + 1); 
                if(b == True){
                    aux_ex = (*(*(p_ex + (i___0_ex - 1)) + k_ex) * *(*(t + k_ex) + j_ex));  
                }
                b = *(m + it);
                it = (it + 1); 
                switch(b){
                    case True:
                        *(*(p_ex + i___0_ex) + j_ex) = aux_ex; 
                    break;
                    case False:
                        
                    break;
                    case Undefined:
                        if((aux_ex > *(*(p_ex + i___0_ex) + j_ex))){
                            *(*(p_ex + i___0_ex) + j_ex) = aux_ex;  
                        }else{
                             
                        }
                    break;
                } 
            } 
        } 
    }
    b = *(m + it);
    it = (it + 1); 
    if(b == True){
        result_ex = 0.000000;  
    }
    for(int j___0_ex = 0; j___0_ex < n;j___0_ex++){
        b = *(m + it);
        it = (it + 1); 
        switch(b){
            case True:
                *(*(p_ex + (n - 1)) + j___0_ex) = result_ex; 
            break;
            case False:
                
            break;
            case Undefined:
                if((*(*(p_ex + (n - 1)) + j___0_ex) < result_ex)){
                    *(*(p_ex + (n - 1)) + j___0_ex) = result_ex;  
                }else{
                     
                }
            break;
        } 
    }
    if(PRINT){
        cout << endl << "Lazy Viterbi" << endl;
        cout << "Result: " << *result << endl;
    }
}

