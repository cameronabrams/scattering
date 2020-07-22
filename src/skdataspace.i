%module libskdataspace
%include carrays.i
%array_functions(double, array);
%array_functions(int, arrayint);
%inline %{
double get_double(double *a, int index) {
	return a[index];
}
%}
%{
#include <stdlib.h>
#include "skdataspace.h"
%}
extern SkDataSpace * NewSkDataSpace ( int n, int nframes, int nkx, int nky, int nkz);
extern void FreeSkDataSpace ( SkDataSpace * sds );
extern double SkDataSpace_setdq ( SkDataSpace * sds, double dq );
extern double * SkDataSpace_getCoordAddr ( SkDataSpace * sds, int dim );
extern double * SkDataSpace_getCellSizeAddr ( SkDataSpace * sds );
extern double * SkDataSpace_getScatLengthAddr ( SkDataSpace * sds );
extern void SkDataSpace_report ( SkDataSpace * sds );
extern int SkDataSpace_setFrame ( SkDataSpace * sds, int f );
extern void SkDataSpace_scaleCoords ( SkDataSpace * sds, int d );
extern void SkDataSpace_updateSk ( SkDataSpace * sds );
extern void SkDataSpace_outputSk ( SkDataSpace * sds );
