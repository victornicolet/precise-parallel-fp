
#ifndef STRUCT_H
#define STRUCT_H

// datatype to return the result of a partial comparison
enum boolean{
    True,
    False,
    Undefined,
    Useless
};

inline boolean uint8toboolean1(uint8_t u){
    if(1 & u){
        if(1<<1 & u){
            True;
        }else{
            False;
        }
    }else{
        if(1<<1 & u){
            Undefined;
        }else{
            Useless;
        }
    }
}

inline boolean uint8toboolean2(uint8_t u){
    if(1<<2 & u){
        if(1<<3 & u){
            True;
        }else{
            False;
        }
    }else{
        if(1<<3 & u){
            Undefined;
        }else{
            Useless;
        }
    }
}

inline boolean uint8toboolean3(uint8_t u){
    if(1<<4 & u){
        if(1<<5 & u){
            True;
        }else{
            False;
        }
    }else{
        if(1<<5 & u){
            Undefined;
        }else{
            Useless;
        }
    }
}

inline boolean uint8toboolean4(uint8_t u){
    if(1<<6 & u){
        if(1<<7 & u){
            True;
        }else{
            False;
        }
    }else{
        if(1<<7 & u){
            Undefined;
        }else{
            Useless;
        }
    }
}

inline uint8_t booleanstouint8_2(boolean b1, boolean b2){
    uint8_t u = 0;
    switch(b1){
        case Useless:
            break;
        case Undefined:
            u = 0x01;
            break;
        case False:
            u = 0x02;
            break;
        case True:
            u = 0x03;
            break;
    }
    switch(b2){
        case Useless:
            break;
        case Undefined:
            u |= 0x04;
            break;
        case False:
            u |= 0x08;
            break;
        case True:
            u |= 0x0c;
            break;
    }
    return u;
}

inline uint8_t booleanstouint8(boolean b1, boolean b2, boolean b3, boolean b4){
    uint8_t u = 0;
    switch(b1){
        case Useless:
            break;
        case Undefined:
            u = 0x01;
            break;
        case False:
            u = 0x02;
            break;
        case True:
            u = 0x03;
            break;
    }
    switch(b2){
        case Useless:
            break;
        case Undefined:
            u |= 0x04;
            break;
        case False:
            u |= 0x08;
            break;
        case True:
            u |= 0x0c;
            break;
    }
    switch(b3){
        case Useless:
            break;
        case Undefined:
            u |= 0x10;
            break;
        case False:
            u |= 0x20;
            break;
        case True:
            u |= 0x30;
            break;
    }
    switch(b4){
        case Useless:
            break;
        case Undefined:
            u |= 0x40;
            break;
        case False:
            u |= 0x80;
            break;
        case True:
            u |= 0xc0;
            break;
    }
    return u;
}

struct memo{
    bool useful1;
    boolean useful2;
};
#endif
