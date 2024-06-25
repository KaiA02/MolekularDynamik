//
// Created by jh on 25.06.2024.
//

#ifndef EPSILONSIGMA_H
#define EPSILONSIGMA_H

class EpsilonSigma{
private:
    int type1;
    int type2;
    double sigma;
    double epsilon;

public:
    inline EpsilonSigma(int t1, int t2, double s, double e){
        type1 = t1;
        type2 = t2;
        sigma = s;
        epsilon = e;
    };
    inline int getT1(){
        return type1;
    };
    inline int getT2(){
        return type2;
    };
    inline double getSigma(){
        return sigma;
    };
    inline double getEpsilon(){
        return epsilon;
    };
    inline bool isRight(int t1, int t2) {
        return (type1 == t1 && type2 == t2) || (type2 == t1 && type1 == t2);
    }

};

#endif //EPSILONSIGMA_H
