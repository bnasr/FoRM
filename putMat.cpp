#include "mat.h"
#include "mex.h"
#include "matrix.h"


#pragma comment(lib, "libmat.lib")
#pragma comment(lib, "libmx.lib")
#pragma comment(lib, "libmex.lib")

extern "C"
{

}


template<class T>
void Array2Mat(T ***Array, char *ArrayName, int nRow, int nCol, char *MatName)
{
    MATFile *pmat;
    
    mxArray *Mat;
    Mat=mxCreateDoubleMatrix(nRow, nCol, mxREAL);

    memcpy(mxGetPr(Mat), *Array, nRow* nCol* sizeof(T));
    
    pmat=matOpen(MatName, "w");
    matPutVariable(pmat, ArrayName, Mat);

    matClose(pmat);
    mxDestroyArray(Mat);
}