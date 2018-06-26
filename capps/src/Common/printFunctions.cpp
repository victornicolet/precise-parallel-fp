#include <iostream>

#include <vector>
#include <emmintrin.h>
#include "printFunctions.hpp"
#include "interval_arithmetic.hpp"

using namespace std;


void printVectorInt(vector<__m128d> v){
    for(vector<__m128d>::iterator it = v.begin(); it != v.end(); it++){
        print(*it);
        cout << endl;
    }
}

void printBoolean(boolean b){
    switch(b){
        case True :
            cout << "True";
            break;
        case False :
            cout << "False";
            break;
        case Undefined :
            cout << "Undefined";
            break;
        case Useless :
            cout << "Useless";
            break;
    }
}    

void printMemo(memo* da, int s){
    for(int i = 0; i < s; i++){
        cout << da[i].useful1 << "," ;
        printBoolean(da[i].useful2);
        cout << endl;
    } 
}
